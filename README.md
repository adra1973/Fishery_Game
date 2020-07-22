# Fishery_Game
Program developed for calculation of a Nash Equilibrium in a Stochastic Dynamic Game
## Description:
This program considers a game among to companies dedicated to capture fish resources in a certain region.

The model includes 11 states of the quantity of fish available in that region.

The two players have four actions to decide the quantity of fish to capture.

The matrices of payoffs for the two players depend of the quantity of fish captured by them.

The matrices of transition consider a Binomial distribution for the disturbances in the quantity of fish in each stage.

Also it is possible run an estimation with the empiric distribution in order to construct the transition matrices.

in order to simplifly the calculations, the actions for both players are the same and thus we have a symmetrical game where the solution is the same for both players.

The Nash Equilibriums are calculated with a procedure using the McKelvey formula to get a function where the minimum points of this function correspond to the Nash Equilibrium.

The results of the program are located in a text file by rows with the number of the stage, the pure strategies and the value of the game for each state.

To run the program, copy the code run it and that is.

You can change some variables and specifications like number of iterations for empiric distribution, number of stages of projection, 
probability of Binomial distribution, discount factor, etc.

I have added comments in every part of the program in order to facilitate its handling.

In case of doubts, suggestions or comments, you can contact its author: Alan Robles at this e-mail adress: alan_daniel@yahoo.com

Thank you.
