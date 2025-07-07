Genint-PWM is a short C program that can generate sequences containing matches to PWMs, or count 
the number of PWM matches in user-given sequences (-file option)

The program is probabilistic and will not give exactly the same results after each run
(it does not use set thresholds, but scales the probabilities described by the PWM and tests 
whether random number is below the set probability)
