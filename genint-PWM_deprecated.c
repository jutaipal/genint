/*  version 0.8                                                                                                                       */
/*  Modified genint-PWM: Can generate sequences OR count matches from input file                                                      */
/*  usage: ./genint-PWM [background PWM] [signal PWM or file with list of PWMs] [number of sequences] [p-value for matches as -log10] [-file input_sequences.txt] */
/*  NEW: Add -file option to count matches in existing sequences instead of generating new ones                                        */
/*  based on original genint-PWM by J. Taipale                                                                                         */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>

#define NUM_STATES 256
static uint64_t xoshiro_states[NUM_STATES][4];
static uint8_t current_state = 0;

#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>

static int has_rdrand(void) {
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        return (ecx & bit_RDRND);
    }
    return 0;
}

static int get_rdrand64(uint64_t* value) {
    unsigned char success;
    __asm__ __volatile__("rdrand %0; setc %1"
                        : "=r" (*value), "=qm" (success)
                        :
                        : "cc");
    return (int)success;
}
#endif

#define NUM_SCALING_ROUNDS 10
static const int target_matches_by_round[NUM_SCALING_ROUNDS] = {3, 10, 33, 100, 333, 1000, 2000, 3000, 4000, 5000};
static const double sampling_error_tolerance = 10.0;
static int pwms_still_scaling;

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t xoshiro256ss(void) {
    uint64_t* state = xoshiro_states[current_state];
    current_state = (current_state + 1) & 0xFF;
    
    const uint64_t result = rotl(state[1] * 5, 7) * 9;
    const uint64_t t = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= t;
    state[3] = rotl(state[3], 45);

    return result;
}

static inline long double xoshiro_to_long_double(void) {
    return (xoshiro256ss() >> 11) * 0x1.0p-53L;
}

short int max_width_of_pwm = 1000;
short int contacts = 0;

struct normalized_pwm {
    char *name;
    char *seed;
    short int width;
    long int max_counts;
    double *information_content;
    short int *original_position;
    double *position_score;
    long int *total_counts_for_column;
    double **fraction;
    short int negative_values_allowed;
    long double scale_factor;
    long matches_in_round;
    long trials_in_round;
    int current_round;
    int scaling_completed;
    double target_frequency;
    long total_final_matches;
};

short int normalized_pwm_init(struct normalized_pwm *i, char *name, short int width, double initial_value) {
    short int maximum_width = max_width_of_pwm;
    short int counter;
    short int counter2;
    (*i).negative_values_allowed = 0;
    (*i).name = malloc(100);
    strcpy((*i).name, name);
    (*i).seed = malloc(1000);
    strcpy((*i).seed, "UNKNOWN");
    (*i).width = width;
    (*i).max_counts = initial_value;
    (*i).fraction = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
    (*i).information_content = malloc(sizeof(double) * maximum_width + 5);
    (*i).position_score = malloc(sizeof(double) * maximum_width + 5);
    (*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
    (*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

    for(counter = 0; counter < 5 + contacts * 12; counter++) {
        (*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
        for(counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
    }
    for(counter2 = 0; counter2 < maximum_width; counter2++) {
        (*i).information_content[counter2] = 0;
        (*i).position_score[counter2] = 0;
        (*i).original_position[counter2] = counter2;
        (*i).total_counts_for_column[counter2] = 0;
    }
    
    (*i).scale_factor = 1.0;
    (*i).matches_in_round = 0;
    (*i).trials_in_round = 0;
    (*i).current_round = 0;
    (*i).scaling_completed = 0;
    (*i).target_frequency = 0.0;
    (*i).total_final_matches = 0.0;
    return(0);
}

short int normalized_pwm_free (struct normalized_pwm *i) {
    short int counter;
    free((*i).name);
    free((*i).information_content);
    free((*i).position_score);
    free((*i).total_counts_for_column);
    for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
    free((*i).fraction);
    return(0);
}

short int Normalize_pwm (struct normalized_pwm *n) {
    short int counter;
    short int position;
    double total_nucleotides = 0;
    double normalized_value = 0;
    for (position = 0; position < (*n).width; position++) {
        for (counter = 0, total_nucleotides = 0; counter < 4; counter++) {
            if ((*n).fraction[counter][position] > 0) total_nucleotides += (*n).fraction[counter][position];
            else if ((*n).negative_values_allowed == 1) total_nucleotides += -(*n).fraction[counter][position];
        }
        
        for (counter = 0; counter < 4; counter++) {
            normalized_value = ((double) (*n).fraction[counter][position]) / total_nucleotides;
            if ((normalized_value < 0) && ((*n).negative_values_allowed == 0)) normalized_value = 0;
            (*n).fraction[counter][position] = normalized_value;
        }
    }
    return (0);
}

short int Load_pwm (struct normalized_pwm *p, char *filename, short int normalize) {
    long int counter;
    char text1;
    short int line = 0;
    short int pwm_position = 0;
    char *current_string;
    current_string = malloc(200);
    FILE *pwmfile;
    if ((pwmfile = fopen(filename, "r")) == (void *)0) {printf("\nNo File: %s", filename); exit (2);}
    for(line = 0; line <= 3 + contacts * 12;) {
        for(counter = 0; counter < 30; counter++) {
            text1 = getc(pwmfile);
            if (text1 == EOF || text1 == '\n' || text1 == '\t') {
                current_string[counter] = '\0';
                if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) {
                    (*p).fraction[line][pwm_position] = atof(current_string);
                    pwm_position++;
                }
                if (text1 == '\n' || text1 == EOF) {(*p).width = pwm_position; line++; pwm_position = 0;}
                break;
            }
            current_string[counter]=text1;
        }
    }
    free (current_string);
    if (normalize == 1) Normalize_pwm(p);
    if (text1 == EOF && line != 3) return(1);
    else return (0);
}

static inline long double Scaled_probability(long double scale) {
    long double p = xoshiro_to_long_double();
    if (p <= 0.0L || p >= 1.0L || scale == 1.0L) {
        return p;
    }
    long double odds = p / (1.0L - p);
    odds *= scale;
    return odds / (1.0L + odds);
}

// RANDOM SEEDER FUNCTIONS
static uint64_t get_system_entropy(uint64_t mixer) {
    uint64_t value;
    static int first_call = 1;  // Only print on first call
    
    // Try /dev/urandom first (works on all Unix systems)
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd != -1) {
        if (read(fd, &value, sizeof(value)) == sizeof(value)) {
            close(fd);
            if (first_call) {
                printf("\tEntropy source: /dev/urandom (hardware RNG)\n");
                first_call = 0;
            }
            return value ^ mixer;
        }
        close(fd);
        if (first_call) {
            printf("\tWarning: /dev/urandom read failed\n");
        }
    } else {
        if (first_call) {
            printf("\tWarning: /dev/urandom open failed\n");
        }
    }
    
    // Fallback 1: Try arc4random (available on BSD/macOS)
    #if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__)
    value = ((uint64_t)arc4random() << 32) | arc4random();
    if (first_call) {
        printf("\tEntropy source: arc4random() (system RNG)\n");
        first_call = 0;
    }
    return value ^ mixer;
    #endif
    
    // Fallback 2: Improved timestamp + process entropy (for other systems)
    if (first_call) {
        printf("\tEntropy source: timestamp fallback (poor quality)\n");
        first_call = 0;
    }
    
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    
    value = ((uint64_t)ts.tv_sec << 32) ^ (uint64_t)ts.tv_nsec;
    value ^= (uint64_t)getpid() << 16;
    value ^= (uint64_t)&value;  // Stack address
    value ^= (uint64_t)clock() << 8;
    
    // Add variable delay based on mixer to break correlation
    usleep((mixer & 0xFF) + 1);  // 1-256 microseconds
    
    // Get another timestamp after delay
    clock_gettime(CLOCK_REALTIME, &ts);
    value ^= ((uint64_t)ts.tv_nsec << 16);
    
    value ^= mixer;
    
    // Simple hash mixing
    value ^= value >> 33;
    value *= 0xff51afd7ed558ccdULL;
    value ^= value >> 33;
    
    return value;
}

void seed_random(void) {
    int use_rdrand = 0;
#if defined(__x86_64__) || defined(__i386__)
    use_rdrand = has_rdrand();
#endif
    
    for (int s = 0; s < NUM_STATES; s++) {
        for (int i = 0; i < 4; i++) {
            uint64_t mixer = ((uint64_t)s << 32) | (i * 0x100000001b3ULL);
            
#if defined(__x86_64__) || defined(__i386__)
            if (use_rdrand) {
                for (int tries = 0; tries < 10; tries++) {
                    if (get_rdrand64(&xoshiro_states[s][i])) {
                        break;
                    }
                }
                if (xoshiro_states[s][i] == 0) {
                    xoshiro_states[s][i] = get_system_entropy(mixer);
                }
            } else
#endif
            {
                xoshiro_states[s][i] = get_system_entropy(mixer);
            }
        }
    }
}

int find_and_improve_best_match(short int* sequence_value, struct normalized_pwm* qs, int num_signal_pwms, int Nlength, int* match_position, long double* final_random_number, long double* final_score, int* matching_pwm, int improve_by_mutation, long int sequences_to_print) {
   long double best_score = 0;
   #define MAX_BEST_POSITIONS 1000
   struct {
       int pos;
       int pwm_index;
       int is_reverse;
   } best_positions[MAX_BEST_POSITIONS];
   int num_best = 0;
   
   for(int p = 0; p < num_signal_pwms; p++) {
       double current_freq = (double)qs[p].total_final_matches / sequences_to_print;
       if(current_freq >= qs[p].target_frequency) continue;

       for (int pos = 0; pos < Nlength-qs[p].width; pos++) {
           long double forward_score = 1;
           for(int i = 0; i < qs[p].width; i++) {
               int nuc = sequence_value[pos+i];
               forward_score *= qs[p].fraction[nuc][i];
           }
           
           if (forward_score >= best_score * 0.999) {
               if (forward_score > best_score * 1.001) {
                   best_score = forward_score;
                   num_best = 0;
               }
               if (num_best < MAX_BEST_POSITIONS) {
                   best_positions[num_best].pos = pos;
                   best_positions[num_best].pwm_index = p;
                   best_positions[num_best].is_reverse = 0;
                   num_best++;
               }
           }
           
           long double reverse_score = 1;
           for(int i = 0; i < qs[p].width; i++) {
               int nuc = sequence_value[pos+i];
               reverse_score *= qs[p].fraction[3-nuc][qs[p].width-1-i];
           }
           
           if (reverse_score >= best_score * 0.999) {
               if (reverse_score > best_score * 1.001) {
                   best_score = reverse_score;
                   num_best = 0;
               }
               if (num_best < MAX_BEST_POSITIONS) {
                   best_positions[num_best].pos = pos;
                   best_positions[num_best].pwm_index = p;
                   best_positions[num_best].is_reverse = 1;
                   num_best++;
               }
           }
       }
   }
   
   if (num_best == 0) return 0;
   
   int chosen = xoshiro256ss() % num_best;
   int best_pos = best_positions[chosen].pos;
   int best_pwm = best_positions[chosen].pwm_index;
   int is_reverse = best_positions[chosen].is_reverse;
   
   if (!improve_by_mutation) {
       long double current_random = Scaled_probability(qs[best_pwm].scale_factor);
       if (current_random < best_score) {
           *match_position = best_pos;
           *final_score = best_score;
           *final_random_number = current_random;
           *matching_pwm = best_pwm;
           return is_reverse ? 2 : 1;
       }
       return 0;
   }
   
    /* ========== ONE-SHOT “WRITE THE PWM” BLOCK — forward & reverse OK ========== */

    {
        /* ---- 1.  draw a fresh motif instance and write it ------------- */
        long double new_score = 1.0L;                          /* product */

        const int width = qs[best_pwm].width;

        for (int rel = 0; rel < width; ++rel) {

            /* Pick the PWM column that corresponds to this sequence index
               -----------------------------------------------------------
               Forward : sequence[best_pos+rel] ↔ PWM column  rel
               Reverse : sequence[best_pos+rel] ↔ PWM column (width-1-rel)
            */
            int pwm_col = is_reverse ? (width - 1 - rel) : rel;

            /* ---- sample a base according to that PWM column ----------- */
            long double r   = xoshiro_to_long_double();
            long double cum = 0.0L;
            int base_fwd    = 0;                    /* 0=A 1=C 2=G 3=T */

            for (; base_fwd < 4; ++base_fwd) {
                cum += qs[best_pwm].fraction[base_fwd][pwm_col];
                if (r < cum) break;
            }
            if (base_fwd == 4) base_fwd = 3;       /* numeric guard    */

            /* ---- write base into sequence ----------------------------- */
            int seq_idx = best_pos + rel;          /* same index both strands */
            sequence_value[seq_idx] = is_reverse ? (3 - base_fwd)   /* complement */
                                                 :  base_fwd;       /* forward    */

            /* ---- accumulate score ------------------------------------- */
            new_score *= qs[best_pwm].fraction[base_fwd][pwm_col];
        }

        /* ---- 2.  acceptance test ------------------------------------- */
        long double thr = Scaled_probability(qs[best_pwm].scale_factor);

        if (thr < new_score) {                     /* motif accepted   */
            *match_position      = best_pos;
            *final_score         = new_score;
            *final_random_number = thr;
            *matching_pwm        = best_pwm;
            return is_reverse ? 2 : 1;             /* 1 = Fwd, 2 = Rev */
        }

        /* ---- 3.  motif rejected – “no match” as original code --------- */
        return 0;
    }
    /* ================== end of replacement block ====================== */
}

int main (int argc, char *argv[]) {
    seed_random();
    
    char *DNA = "ACGTN";
    char *dnalc = "acgtn";
    short int Nlength;
    int pwms_still_scaling = 0;
    int nomutate = 0;
    
    // NEW: Check for -file option
    char *input_sequence_file = NULL;
    int use_input_file = 0;
    
    // Check for -file flag
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-file") == 0 && i + 1 < argc) {
            input_sequence_file = argv[i + 1];
            use_input_file = 1;
            // Shift remaining arguments left
            for(int j = i; j < argc - 2; j++) {
                argv[j] = argv[j + 2];
            }
            argc -= 2;
            break;
        }
    }
    
    // Check for -nomutate flag
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-nomutate") == 0) {
            nomutate = 1;
            for(int j = i; j < argc-1; j++) {
                argv[j] = argv[j+1];
            }
            argc--;
            break;
        }
    }
    
    if(argc != 5) {
        fprintf(stderr, "\ngenint-PWM v0.23 with file input option\n\n");
        fprintf(stderr, "Usage: genint-PWM [-nomutate] <background PWM> <signal PWM/PWM list> <number of sequences> <-log10 p-value> [-file input_sequences.txt]\n\n");
        fprintf(stderr, "Generate sequences with PWM matches OR count matches in existing sequences\n");
        fprintf(stderr, "  -nomutate             : Use pure rejection sampling without mutation improvement\n");
        fprintf(stderr, "  -file <sequences.txt> : Count matches in existing sequences instead of generating\n");
        fprintf(stderr, "  <background PWM>      : PWM file for background nucleotide frequencies\n");
        fprintf(stderr, "  <signal PWM/PWM list> : Single PWM file or text file with list of PWM files and percentages\n");
        fprintf(stderr, "  <number of sequences> : Number of sequences to generate/analyze\n");
        fprintf(stderr, "  <-log10 p-value>      : Desired p-value as -log10 (e.g., 5 for p=10^-5)\n\n");
        if(use_input_file) {
            fprintf(stderr, "FILE MODE: Will count matches in sequences from %s\n", input_sequence_file);
        } else {
            fprintf(stderr, "GENERATE MODE: Will generate new sequences\n");
        }
        exit(1);
    }
    
    char *backgroundPWM_name = malloc(1000);
    strcpy(backgroundPWM_name, argv[1]);
    
    char *signalPWM_name = malloc(1000);
    strcpy(signalPWM_name, argv[2]);
    
    long int number_of_generated_sequences = atoi(argv[3]);
    long int sequences_to_print = number_of_generated_sequences;
    
    int minus_log10_pvalue = atoi(argv[4]);
    long double target_pvalue = pow(10, -minus_log10_pvalue);
    
    if(use_input_file) {
        printf("\tFILE MODE: Counting matches in %s\n", input_sequence_file);
    } else {
        printf("\tGENERATE MODE: Generating sequences\n");
    }
    
    // Load background PWM
    struct normalized_pwm qb;
    normalized_pwm_init(&qb, "empty", Nlength * 2, 0);
    Load_pwm (&qb, backgroundPWM_name, 1);
    strcpy(qb.name, backgroundPWM_name);
    Normalize_pwm(&qb);
    Nlength = qb.width+1;
    
    // Load signal PWMs
    struct normalized_pwm *qs;
    int num_signal_pwms = 1;
    
    // Check if file is PWM or list
    FILE *test_file = fopen(signalPWM_name, "r");
    char test_line[1000];
    if (fgets(test_line, sizeof(test_line), test_file)) {
        for(int i = 0; test_line[i] != '\0'; i++) {
            if (isalpha(test_line[i])) {
                num_signal_pwms = 1;
                while(fgets(test_line, sizeof(test_line), test_file)) {
                    num_signal_pwms++;
                }
                break;
            }
        }
    }
    fclose(test_file);
    
    qs = malloc(num_signal_pwms * sizeof(struct normalized_pwm));
    
    if (num_signal_pwms == 1) {
        normalized_pwm_init(&qs[0], "empty", Nlength * 2, 0);
        Load_pwm(&qs[0], signalPWM_name, 1);
        strcpy(qs[0].name, signalPWM_name);
        Normalize_pwm(&qs[0]);
        qs[0].target_frequency = 0.05;
    } else {
        FILE *list_file = fopen(signalPWM_name, "r");
        char line[1000];
        int pwm_index = 0;
        while(fgets(line, sizeof(line), list_file)) {
            char pwm_filename[1000];
            double percentage;
            if(sscanf(line, "%s %lf", pwm_filename, &percentage) == 2) {
                normalized_pwm_init(&qs[pwm_index], "empty", Nlength * 2, 0);
                Load_pwm(&qs[pwm_index], pwm_filename, 1);
                strcpy(qs[pwm_index].name, pwm_filename);
                Normalize_pwm(&qs[pwm_index]);
                qs[pwm_index].target_frequency = percentage / 100.0;
                pwm_index++;
            }
        }
        fclose(list_file);
    }
    
    // Set scale factors using original logic but simplified
    for(int i = 0; i < num_signal_pwms; i++) {
        long double log_expected_score = 0.0;
        
        // Calculate expected score under uniform background (0.25 each base)
        for(int pos = 0; pos < qs[i].width; pos++) {
            long double pos_score = 0.0;
            for(int base = 0; base < 4; base++) {
                pos_score += 0.25 * qs[i].fraction[base][pos];  // Uniform background
            }
            if(pos_score > 0) {
                log_expected_score += log(pos_score);
            }
        }
        
        // Original genint-PWM scale factor calculation
        long double log_target = log(target_pvalue);
        qs[i].scale_factor = exp(log_expected_score - log_target);
        
        printf("\tPWM %s: expected_score=%.2e, scale_factor=%.6Lf\n",
               qs[i].name, (double)exp(log_expected_score), qs[i].scale_factor);
    }
    
    long int current_sequence_position;
    long double current_random_number;
    double cutoff;
    short int base;
    short int nucleotide;
    short int pwm_position;
    long int round = 0;
    
    long double score;
    short int win = 0;
    short int pos;
    long int total_matches = 0;
    
    short int *sequence_value = malloc(Nlength * sizeof(short int));
    if (!sequence_value) { perror("malloc"); exit(3); }
    long int total_number_of_test_sequences_generated_for_printing = 0;
    short int add_random_sequence;
    
    if(use_input_file) {
        // FILE MODE: Read sequences from input file
        FILE *seq_file = fopen(input_sequence_file, "r");
        if(!seq_file) {
            printf("Error: Cannot open sequence file %s\n", input_sequence_file);
            return 1;
        }
        
        char line[2048];
        int *match_counts = calloc(num_signal_pwms, sizeof(int));
        int total_sequences_processed = 0;
        
        printf("\nProcessing sequences from file...\n");
        
        while(fgets(line, sizeof(line), seq_file) && number_of_generated_sequences > 0) {
            // Clean up line
            int len = strlen(line);
            while(len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ' || line[len-1] == '\t')) {
                line[len-1] = '\0';
                len--;
            }
            
            // Skip empty lines, FASTA headers, and tab-separated lines
            if(len == 0 || line[0] == '>' || strchr(line, '\t') != NULL) continue;
            
            // Validate DNA sequence
            int valid_sequence = 1;
            for(int i = 0; i < len; i++) {
                char c = toupper(line[i]);
                if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                    valid_sequence = 0;
                    break;
                }
            }
            if(!valid_sequence) continue;
            
            // Convert sequence to numeric array
            int actual_len = (len < max_width_of_pwm) ? len : max_width_of_pwm - 1;
            for(int i = 0; i < actual_len; i++) {
                char c = toupper(line[i]);
                switch(c) {
                    case 'A': sequence_value[i] = 0; break;
                    case 'C': sequence_value[i] = 1; break;
                    case 'G': sequence_value[i] = 2; break;
                    case 'T': sequence_value[i] = 3; break;
                    default: sequence_value[i] = 0; break;
                }
            }
            
            total_sequences_processed++;
            
            // Test each PWM against this sequence
            for(int pwm_idx = 0; pwm_idx < num_signal_pwms; pwm_idx++) {
                if(actual_len < qs[pwm_idx].width) continue;
                
                int found_match = 0;
                
                // Check forward and reverse orientations
                for(int pos = 0; pos <= actual_len - qs[pwm_idx].width && !found_match; pos++) {
                    // Forward orientation
                    long double forward_score = 1.0;
                    for(int i = 0; i < qs[pwm_idx].width; i++) {
                        forward_score *= qs[pwm_idx].fraction[sequence_value[pos + i]][i];
                    }
                    
                    long double threshold = Scaled_probability(qs[pwm_idx].scale_factor);
                    if(threshold < forward_score) {
                        found_match = 1;
                        break;
                    }
                    
                    // Reverse complement
                    long double reverse_score = 1.0;
                    for(int i = 0; i < qs[pwm_idx].width; i++) {
                        int complement = 3 - sequence_value[pos + i];
                        reverse_score *= qs[pwm_idx].fraction[complement][qs[pwm_idx].width - 1 - i];
                    }
                    
                    threshold = Scaled_probability(qs[pwm_idx].scale_factor);
                    if(threshold < reverse_score) {
                        found_match = 1;
                        break;
                    }
                }
                
                if(found_match) {
                    match_counts[pwm_idx]++;
                }
            }
            
            number_of_generated_sequences--;
            
            //if(total_sequences_processed % 1000 == 0) {
            //    printf("Processed %d sequences...\n", total_sequences_processed);
            //}
        }
        
        fclose(seq_file);
        
        // Output results in format suitable for genint-PWM input
        printf("\nMatch results:\n");
        FILE *results_file = fopen("match_results.txt", "w");
        for(int i = 0; i < num_signal_pwms; i++) {
            double percentage = (total_sequences_processed > 0) ? (100.0 * match_counts[i] / total_sequences_processed) : 0.0;
            printf("PWM %s: %d/%d matches (%.2f%%)\n", qs[i].name, match_counts[i], total_sequences_processed, percentage);
            if(results_file) {
                fprintf(results_file, "%s\t%.2f\n", qs[i].name, percentage);
            }
        }
        if(results_file) {
            fclose(results_file);
            printf("Results written to match_results.txt\n");
        }
        
        free(match_counts);
        
    } else {
        // ORIGINAL GENERATION MODE - unchanged
        for (; number_of_generated_sequences > 0;) {
            round++;
            for (current_sequence_position = 0; current_sequence_position < Nlength-1; current_sequence_position++) {
                current_random_number = xoshiro_to_long_double();
                for(cutoff = 0, base = 0; base < 4; base++) {
                    cutoff += qb.fraction[base][current_sequence_position];
                    if (current_random_number < cutoff) break;
                }
                sequence_value[current_sequence_position] = base;
            }
            
            total_number_of_test_sequences_generated_for_printing++;
            
            int match_position;
            add_random_sequence = 0;
            int matching_pwm = -1;
            
            win = find_and_improve_best_match(sequence_value, qs, num_signal_pwms,
                                            Nlength, &match_position, &current_random_number,
                                            &score, &matching_pwm, (nomutate ? 0 : 1), sequences_to_print);
            if (win != 0) {
                pos = match_position;
            }
            
            if(win != 0) {
                double current_freq = (double)qs[matching_pwm].total_final_matches / sequences_to_print;
                if(current_freq >= qs[matching_pwm].target_frequency) {
                    win = 0;
                }
                else {
                    qs[matching_pwm].total_final_matches++;
                    total_matches++;
                }
            }
            
            if (win == 0) {
                double total_target = 0;
                for(int i = 0; i < num_signal_pwms; i++) {
                    total_target += qs[i].target_frequency;
                }
                
                double current_match_seq_freq = (double)(total_matches + 1) / sequences_to_print;
                
                if(current_match_seq_freq > total_target) {
                    add_random_sequence = 1;
                }
            }
            
            // Print sequences (both matching and random)
            if (win != 0 || add_random_sequence == 1) {
                printf("\n");
                for (current_sequence_position = 0; current_sequence_position < Nlength-1; current_sequence_position++) {
                    nucleotide = sequence_value[current_sequence_position];
                    printf("%c",DNA[nucleotide]);
                }
                if (win != 0) {
                    printf("\t%s\t%i%c\t%li\t%.10Lf\t%.10Lf",
                           qs[matching_pwm].name, pos, "FR"[win-1], round,
                           current_random_number, score);
                } else {
                    printf("\tRandom non-matching sequence");
                }
                number_of_generated_sequences--;
            }
        }

        printf("\n\tFinal match statistics:\n");
        printf("\tPWM\tTarget%%\tObserved%%\tMatches\n");
        printf("\t----------------------------------------\n");
        for(int i = 0; i < num_signal_pwms; i++) {
            double observed_percent = 100.0 * qs[i].total_final_matches / sequences_to_print;
            printf("\t%s\t%.1f\t%.1f\t%ld\n",
                   qs[i].name,
                   qs[i].target_frequency * 100.0,
                   observed_percent,
                   qs[i].total_final_matches);
        }
        printf("\t----------------------------------------\n");
        printf("\tTotal sequences with matches: %ld\n", total_matches);
        printf("\tTotal sequences tried: %ld\n", total_number_of_test_sequences_generated_for_printing);
    } // End of else block for generation mode
    
    printf("\n");
    
    free(backgroundPWM_name);
    free(signalPWM_name);
    free(qs);
    return 0;
} // End of main function
