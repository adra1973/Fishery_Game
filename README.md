# Fishery_Game
Program developed for calculation of a Nash Equilibrium in a Stochastic Dynamic Game

## Remark:

This program is used for the numerical calculation in Section 4 of the paper:

[1] Alan D. Robles-Aguilar, David González-Sánchez and J. Adolfo Minjárez-Sosa (2020) *Empirical approximation of Nash equilibria in finite Markov games with discounted payoffs*.(Submited for publication).

### Description:

This program provides a computational algorithm to approximate Nash equilibria for a stochastic finite version of the Great Fish War game, in the settings of paper [1]. In perticular the random disturbances are assumed with Binomial distribution.

In order to simplify the calculations, the action sets for both players are equal. Thus, we have a symmetrical game where the equilibrium strategy for players coincides.

The Nash equilibria are calculated with a procedure that involves the McKelvey formula and the empirical distribution as an estimator of the Binomial distribution.

Specifically, the program provides symmetric equilibrium strategies in the full-information game, as well as symmetric equilibrium strategies for some simulated empirical games, both in the set of mixed strategies.

Among the variables that must be specified for the operation of the program are the number of iterations, the parameters of the Binomial distribution, the discount factor, etc. In fact, to facilitate his handling, we have added comments in each part of the program.

For more information, suggestions or comments contact to Alan Robles at e-mail: alan_daniel@yahoo.com

Thank you
