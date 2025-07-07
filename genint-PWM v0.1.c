/*  generative part of interpretable model v0.1: script that generates sequences from PWMs by rejection sampling                                     */
/*  usage: ./genint-PWM [background PWM] [signal PWM or file with list of PWMs] [number of sequences generated] [p-value for matches as -log10]      */
/*  list must contain text file in which each row contains PWM1filename and percent of matching sequences generated for that PWM, separated by tab   */
/*  written Tuesday 24 Dec 2024 by J. Taipale  (with help from Claude and ChatGPT)                                                                   */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>

// Global constants for scaling
#define NUM_SCALING_ROUNDS 10
static const int target_matches_by_round[NUM_SCALING_ROUNDS] = {3, 10, 33, 100, 333, 1000, 2000, 3000, 4000, 5000};
static const double sampling_error_tolerance = 10.0;  // 10%
static int pwms_still_scaling;

static uint64_t xoshiro_state[4];

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t xoshiro256ss(void) {
    const uint64_t result = rotl(xoshiro_state[1] * 5, 7) * 9;
    const uint64_t t = xoshiro_state[1] << 17;

    xoshiro_state[2] ^= xoshiro_state[0];
    xoshiro_state[3] ^= xoshiro_state[1];
    xoshiro_state[1] ^= xoshiro_state[2];
    xoshiro_state[0] ^= xoshiro_state[3];

    xoshiro_state[2] ^= t;
    xoshiro_state[3] = rotl(xoshiro_state[3], 45);

    return result;
}

static inline long double xoshiro_to_long_double(void) {
    return (xoshiro256ss() >> 11) * 0x1.0p-53L;
}

/* NORMALIZED PWM STRUCTURE */
short int max_width_of_pwm = 1000;
short int contacts = 0;

struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; short int negative_values_allowed; long double scale_factor; long matches_in_round; long trials_in_round; int current_round; int scaling_completed; double target_frequency; long total_final_matches;};
short int normalized_pwm_init(struct normalized_pwm *i, char *name, short int width, double initial_value)
{
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

    for(counter = 0; counter < 5 + contacts * 12; counter++)
    {
        (*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
        for(counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
    }
    for(counter2 = 0; counter2 < maximum_width; counter2++)
    {
        (*i).information_content[counter2] = 0;
        (*i).position_score[counter2] = 0;
        (*i).original_position[counter2] = counter2;
        (*i).total_counts_for_column[counter2] = 0;
    }
    
    // Initialize scaling fields
    (*i).scale_factor = 1.0;
    (*i).matches_in_round = 0;
    (*i).trials_in_round = 0;
    (*i).current_round = 0;
    (*i).scaling_completed = 0;
    (*i).target_frequency = 0.0;
    (*i).total_final_matches = 0.0;
    return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}


/* SUBROUTINE THAT RENORMALIZES NORMALIZED PWM (ROWS IN EACH COLUMN ADD TO 1) */
short int Normalize_pwm (struct normalized_pwm *n)
{
    short int counter;
    short int position;
    double total_nucleotides = 0;
    double normalized_value = 0;
    // printf("\nWidth: %d", (*n).width);
    for (position = 0; position < (*n).width; position++)
    {
        for (counter = 0, total_nucleotides = 0; counter < 4; counter++)
        {
        if ((*n).fraction[counter][position] > 0) total_nucleotides += (*n).fraction[counter][position];
        else if ((*n).negative_values_allowed == 1) total_nucleotides += -(*n).fraction[counter][position];
        }
        
        for (counter = 0; counter < 4; counter++)
        {
        normalized_value = ((double) (*n).fraction[counter][position]) / total_nucleotides;
        if ((normalized_value < 0) && ((*n).negative_values_allowed == 0)) normalized_value = 0;
        (*n).fraction[counter][position] = normalized_value;
        }
    }
    return (0);
}


/* SUBROUTINE THAT LOADS A PWM AND NORMALIZES IT */
short int Load_pwm (struct normalized_pwm *p, char *filename, short int normalize)
{
long int counter;
char text1;
short int line = 0;
short int pwm_position = 0;
char *current_string;
current_string = malloc(200);
FILE *pwmfile;
if ((pwmfile = fopen(filename, "r")) == (void *)0) {printf("\nNo File: %s", filename); exit (2);}
for(line = 0; line <= 3 + contacts * 12;)
{
    for(counter = 0; counter < 30; counter++)
    {
        text1 = getc(pwmfile);
        if (text1 == EOF || text1 == '\n' || text1 == '\t')
        {
        current_string[counter] = '\0';
        if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-'))
        {
        (*p).fraction[line][pwm_position] = atof(current_string);
        //printf("\n%f", (*p).fraction[line][pwm_position]);
        pwm_position++;
        }
        if (text1 == '\n' || text1 == EOF) {(*p).width = pwm_position; line++; pwm_position = 0;}
        break;
        }
        current_string[counter]=text1;
        /* printf ("%c", text1); */
    }
}
free (current_string);
if (normalize == 1) Normalize_pwm(p);
// for (line = 0; line < 4 + contacts * 12; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
if (text1 == EOF && line != 3) return(1);
else return (0);
}

static inline long double Scaled_probability(long double scale)
{
    long double p = xoshiro_to_long_double();  // in (0,1)
    if (p <= 0.0L || p >= 1.0L || scale == 1.0L) {
        return p;
    }
    long double odds = p / (1.0L - p);
    odds *= scale;
    return odds / (1.0L + odds);
}

int check_PWM_match(short int* sequence_value, struct normalized_pwm* qs, int Nlength, int* match_position, long double *final_random_number, long double* final_score) {
    int current_sequence_position;
    int pwm_position;
    int nucleotide;
    long double score;
    long double current_random_number;
    
    // Check forward orientation
    for (current_sequence_position = 0; current_sequence_position < Nlength-qs->width; current_sequence_position++) {
        score = 1;
        long double current_random_number = Scaled_probability(qs->scale_factor);
        for(pwm_position = 0; pwm_position < qs->width && current_random_number < score; pwm_position++) {
            int nucleotide = sequence_value[current_sequence_position+pwm_position];
            score *= qs->fraction[nucleotide][pwm_position];
        }
        
            if (current_random_number < score) {
                *match_position = current_sequence_position;
                *final_score = score;  // Store the winning score
                *final_random_number = current_random_number;
                return 1;  // 1 indicates forward match
            }
    }

    // Check reverse complement
    for (current_sequence_position = 0; current_sequence_position < Nlength-qs->width; current_sequence_position++) {
        score = 1;
        long double current_random_number = Scaled_probability(qs->scale_factor);
        for(pwm_position = 0; pwm_position < qs->width && current_random_number < score; pwm_position++) {
            nucleotide = sequence_value[current_sequence_position+pwm_position];
            score *= qs->fraction[3-nucleotide][qs->width-1-pwm_position];
        }
        
            if (current_random_number < score) {
                *match_position = current_sequence_position;
                *final_score = score;  // Store the winning score
                //printf("\nReverse match!!");
                return 2;  // 2 indicates reverse match
        }
    }
    
    *match_position = -1;
    return 0;  // 0 indicates no match
}

void set_PWM_scales(struct normalized_pwm* qs, int num_pwms, struct normalized_pwm* qb, long double target_pvalue) {
    for(int i = 0; i < num_pwms; i++) {
        long double log_expected_score = 0.0;  // Using sum of logs instead of product
        
        // For each position in PWM
        for(int pos = 0; pos < qs[i].width; pos++) {
            long double pos_score = 0.0;
            for(int base = 0; base < 4; base++) {
                pos_score += qb->fraction[base][0] * qs[i].fraction[base][pos];
            }
            log_expected_score += log(pos_score);  // Sum logs instead of multiplying
        }
        
        // Convert target p-value to log space too
        long double log_target = log(target_pvalue);
        
        // Scale factor in log space
        qs[i].scale_factor = exp(log_expected_score - log_target);
        
        printf("\tPWM: %s\n", qs[i].name);
        printf("\tLog expected score: %.10Lf\n", log_expected_score);
        printf("\tScale factor: %.10Lf\n", qs[i].scale_factor);
        printf("\tTarget: %.10Lf\n\n", target_pvalue);
    }
}

double adjust_PWM_scaling(struct normalized_pwm *pwm) {
    if(pwm->scaling_completed) return 0.0;
    
    double observed_freq = (double)pwm->matches_in_round / pwm->trials_in_round;
    double se = sqrt((observed_freq * (1-observed_freq)) / pwm->trials_in_round);
    double relative_error = (se / observed_freq) * 100.0;
    
    // If we have enough matches for this round
    if(pwm->matches_in_round >= target_matches_by_round[pwm->current_round]) {
        double ratio = pwm->target_frequency / observed_freq;
        
        // Adjust scale factor
        pwm->scale_factor *= ratio;
        
        // Check if final round and converged
        if(pwm->current_round == NUM_SCALING_ROUNDS - 1) {
            if(relative_error < sampling_error_tolerance &&
               fabs(ratio - 1.0) < sampling_error_tolerance/100.0) {
                pwm->scaling_completed = 1;
            }
        } else {
            // Move to next round
            pwm->current_round++;
        }
        
        // Reset counters for next round
        pwm->matches_in_round = 0;
        pwm->trials_in_round = 0;
        
        printf("PWM %s Round %d: matches=%ld trials=%ld error=%.2f%% ratio=%.3f new_scale=%.4Lf\n",
               pwm->name, pwm->current_round, pwm->matches_in_round,
               pwm->trials_in_round, relative_error, ratio, pwm->scale_factor);
    }
    
    return relative_error;
}


/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main (int argc, char *argv[])

{
    
    srand(time(NULL));
    for(int i = 0; i < 4; i++) {
        xoshiro_state[i] = ((uint64_t)rand() << 32) | rand();
    }
    
    char *DNA = "ACGTN";
    char *dnalc = "acgtn";
    short int Nlength;
    
    int pwms_still_scaling;  // Will be initialized after loading PWMs
    
    if(argc != 5 && argc != 6) {
        fprintf(stderr, "\ngenint-PWM v0.21\n\n");
        fprintf(stderr, "Usage: genint-PWM <background PWM> <signal PWM/PWM list> <number of sequences> <-log10 p-value> <optional: %% of matches for single signal PWM>\n\n");
        fprintf(stderr, "Generate sequences with PWM matches at specified p-value\n");
        fprintf(stderr, "  <background PWM>     : PWM file for background nucleotide frequencies\n");
        fprintf(stderr, "  <signal PWM/PWM list>: Single PWM file or text file with list of PWM files and percent of matches for that PWM (no %% sign, tab separated)\n");
        fprintf(stderr, "  <number of sequences>: Number of sequences to generate\n");
        fprintf(stderr, "  <-log10 p-value>     : Desired p-value as -log10 (e.g., 5 for p=10^-5)\n");
        fprintf(stderr, "  <%% of matches>       : Optional: Desired percentage of matches for single signal PWM file (defaults to 5)\n\n");
        exit(1);
    }
    
    char *backgroundPWM_name;
    backgroundPWM_name = malloc(1000);
    strcpy(backgroundPWM_name, argv[1]);
    
    char *signalPWM_name;
    signalPWM_name = malloc(1000);
    strcpy(signalPWM_name, argv[2]);
    
    long int number_of_generated_sequences;
    number_of_generated_sequences = atoi(argv[3]);
    long int sequences_to_print = number_of_generated_sequences;
    
    int minus_log10_pvalue = atoi(argv[4]);  // e.g., 5 for 10^-5
    long double target_pvalue = pow(10, -minus_log10_pvalue);
    
    /* BACKGROUND PWM STRUCTURE */
    struct normalized_pwm qb;
    normalized_pwm_init(&qb, "empty", Nlength * 2, 0);
    
    Load_pwm (&qb, backgroundPWM_name, 1);
    strcpy(qb.name, backgroundPWM_name);
    Normalize_pwm(&qb);
    Nlength = qb.width+1;
    
    /* SIGNAL PWM STRUCTURE */
    struct normalized_pwm *qs;  // Changed to pointer for array
    int num_signal_pwms = 1;   // Default for single PWM case
    
    // Check if file is PWM or list
    FILE *test_file = fopen(signalPWM_name, "r");
    char test_line[1000];
    if (fgets(test_line, sizeof(test_line), test_file)) {
        for(int i = 0; test_line[i] != '\0'; i++) {
            if (isalpha(test_line[i])) {
                // Count number of PWMs in list
                num_signal_pwms = 1;  // Count first line
                while(fgets(test_line, sizeof(test_line), test_file)) {
                    num_signal_pwms++;
                }
                break;
            }
        }
    }
    fclose(test_file);
    
    // Allocate array of PWMs
    qs = malloc(num_signal_pwms * sizeof(struct normalized_pwm));
    
    if (num_signal_pwms == 1) {
        // Single PWM file
        normalized_pwm_init(&qs[0], "empty", Nlength * 2, 0);
        Load_pwm(&qs[0], signalPWM_name, 1);
        strcpy(qs[0].name, signalPWM_name);
        Normalize_pwm(&qs[0]);
        if (argc == 6) qs[0].target_frequency = atof(argv[5]);
        else qs[0].target_frequency = 0.05;  // Default 5% if no percentage given
    } else {
        // List of PWM files with percentages
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
    
    //long double target_pvalue = 1e-5;  // or from argv if you prefer
    set_PWM_scales(qs, num_signal_pwms, &qb, target_pvalue);
    
    //pwms_still_scaling = num_signal_pwms;
    pwms_still_scaling = 0;
    
    signed short int use_background_position = 0;
    long int current_sequence_position;
    long double current_random_number;
    double cutoff;
    short int base;
    short int nucleotide;
    short int revcomp;
    short int pwm_position;
    long int round = 0;
    
    long double score;
    short int win = 0;
    short int pos;
    long int total_matches = 0;
    
    short int sequence_value[Nlength+1];
    long double old_scale;
    double observed_freq;
    double ratio;
    long int total_number_of_test_sequences_generated_for_printing = 0;
    short int add_random_sequence;
    
    for (; number_of_generated_sequences > 0;)
    {
        round++;
        for (current_sequence_position = 0; current_sequence_position < Nlength-1; current_sequence_position++)
        {
            //current_random_number = (long double) rand () / (long double) RAND_MAX;
            current_random_number = xoshiro_to_long_double();
            for(cutoff = 0, base = 0; base < 4; base++)
            {
                cutoff += qb.fraction[base][current_sequence_position];
                if (current_random_number < cutoff) break;
            }
            sequence_value[current_sequence_position] = base;
        }
        
        total_number_of_test_sequences_generated_for_printing++; // will be zeroed when scaling is completed for all PWMs
        
        int match_position;
        add_random_sequence = 0;
        int matching_pwm = -1;  // To store which PWM matched
        for(int i = 0; i < num_signal_pwms; i++) {
            win = check_PWM_match(sequence_value, &qs[i], Nlength, &match_position, &current_random_number, &score);
            if (win != 0) {
                pos = match_position;
                matching_pwm = i;
                break;
            }
        }
        
        // Count trial for all PWMs still scaling
        if(pwms_still_scaling > 0) {
            for(int i = 0; i < num_signal_pwms; i++) {
                if(!qs[i].scaling_completed) {
                    qs[i].trials_in_round++;
                }
            }
        }
        
        if(win != 0) {
            if(pwms_still_scaling > 0) {  // Still in scaling phase
                
                if(!qs[matching_pwm].scaling_completed) {
                    qs[matching_pwm].matches_in_round++;
                    
                    if(qs[matching_pwm].matches_in_round >= target_matches_by_round[qs[matching_pwm].current_round]) {
                        observed_freq = (double) qs[matching_pwm].matches_in_round / (double)
                        qs[matching_pwm].trials_in_round;
                        
                        ratio = observed_freq / qs[matching_pwm].target_frequency;
                        
                        if(ratio > 100.0) ratio = 100.0;
                        if(ratio < 0.01) ratio = 0.01;
                        
                        old_scale = qs[matching_pwm].scale_factor;
                        qs[matching_pwm].scale_factor *= ratio;
                        
                        // if(qs[matching_pwm].scale_factor < 1e-9) { qs[matching_pwm].scale_factor = 1e-9;}
                        
                        printf("PWM %s round %d: matches %ld trials %ld freq %.4f expected %.4f adjustment %.4f scale %.12Lf -> %.12Lf\n",qs[matching_pwm].name, qs[matching_pwm].current_round,qs[matching_pwm].matches_in_round, qs[matching_pwm].trials_in_round,observed_freq * 100, qs[matching_pwm].target_frequency * 100, ratio, old_scale, qs[matching_pwm].scale_factor);
                        
                        if(qs[matching_pwm].current_round == NUM_SCALING_ROUNDS - 1) {
                            qs[matching_pwm].scaling_completed = 1;
                            pwms_still_scaling--;
                            if (pwms_still_scaling == 0) total_number_of_test_sequences_generated_for_printing = 0;
                            printf("PWM %s scaling completed, observed freq %.2f%% vs %.2f%% expected\n", qs[matching_pwm].name, observed_freq * 100, qs[matching_pwm].target_frequency * 100);
                        } else {
                            qs[matching_pwm].current_round++;
                        }
                        
                        qs[matching_pwm].matches_in_round = 0;
                        qs[matching_pwm].trials_in_round = 0;
                    }
                }
            }
            
        }

            //printf("After quota check: win=%d add_random=%d\n", win, add_random_sequence);

            // Check if we need random sequence (when no PWM match)
            if (win == 0) {
                double total_target = 0;
                for(int i = 0; i < num_signal_pwms; i++) {
                    total_target += qs[i].target_frequency;
                }
                
                // Calculate based on target sequences and actual accepted matches
                double current_match_seq_freq = (double)(total_matches + 1) / sequences_to_print;
                
                /*printf("Random check: accepted_matches=%ld target_seqs=%ld curr_match_freq=%.4f max_match_target=%.4f\n",
                       total_matches, number_of_generated_sequences,
                       current_random_freq, total_target); */
                
                if(current_match_seq_freq > total_target) {
                    add_random_sequence = 1;
                }
            }



            if (win != 0) {
                // Then check if PWM match is over quota
                /* printf("Before quota check: win=%d PWM=%d matches=%ld seqs=%ld\n",
                      win, matching_pwm, qs[matching_pwm].total_final_matches, number_of_generated_sequences); */
               double current_freq = (double)qs[matching_pwm].total_final_matches / sequences_to_print;
               if(current_freq >= qs[matching_pwm].target_frequency) {
                   win = 0;  // Already met quota for this PWM
               }
                   else {
                       qs[matching_pwm].total_final_matches++; // counts matches
                       total_matches++;
                   }

            }
            
            /* printf("Final decision: win=%d add_random=%d\n", win, add_random_sequence); */
            
            
            // Print phase
            if (pwms_still_scaling == 0 && (win != 0 || add_random_sequence == 1))
            {
                printf("\n");
                for (current_sequence_position = 0; current_sequence_position < Nlength-1; current_sequence_position++) {
                    nucleotide = sequence_value[current_sequence_position];
                    printf("%c",DNA[nucleotide]);
                }
                if (win != 0) printf("\t%s\t%i%c\t%li\t%.10Lf\t%.10Lf", qs[matching_pwm].name, pos, "FR"[win-1], round, current_random_number, score);
                else printf ("\tRandom non-matching sequence");
                number_of_generated_sequences--;
            }


        } // End of sequence generation loop

        printf("\n\tFinal match statistics:\n");
        printf("\tPWM\tTarget%%\tObserved%%\tMatches\n");
        printf("\t----------------------------------------\n");
    for(int i = 0; i < num_signal_pwms; i++) {
        // Change this line:
        double observed_percent = 100.0 * qs[i].total_final_matches /
                                sequences_to_print;  // Use initial target number
        printf("\t%s\t%.1f\t%.1f\t%ld\n",
               qs[i].name,
               qs[i].target_frequency * 100.0,
               observed_percent,
               qs[i].total_final_matches);
    }
        printf("\t----------------------------------------\n");
        printf("\tTotal sequences with matches: %ld\n", total_matches);
        printf("\tTotal sequences tried: %ld\n", total_number_of_test_sequences_generated_for_printing);
        printf("\n");
        }



