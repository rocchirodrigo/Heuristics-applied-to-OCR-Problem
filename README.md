# Conclusion Thesis Work
Tabu Search Algorithm is developted in order to find a set of optimal values for Plug Settings (PS) and Dial Setting for a given number 
of Overcurrent Relays (OCRs). The performance was compared to simpler methods as Linear Programming (given PS lowest possible values, the
optimal Dial setting is obtained) and Heuristic (sort of a Hill Climbing algorithm, finds optimal values for both PS and Dial).

Tabu Search is an algorithm which can be compared to the Hill Climbing, with the exception that once an optimal value is found, the program 
is prohibited to search in that point again, being compelled to explore the Search Space to avoid local optima.

Further improvements would be utilizing Diversification and Intensification methods combined to find better values or to reduce execution 
time of the algorithm.
