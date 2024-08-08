# STRING MAPPING (subsequence alignment problem)

### Goal: The goal of this assignment is to take a complex new problem and formulate and solve it as search. Formulation as search is an integral skill of AI that will come in handy whenever you are faced with a new problem. Heuristic search will allow you to find optimal solutions. Local search may not find the optimal solution, but is usually able to find good solutions for really large problems.

#### Scenario: You are a Genetics researcher working on phylogeny data – you have gene sequences of various organisms. You want to prove that some organisms are more related than others. You decide to map the strings onto each other and calculate mapping scores between each pair of organisms. The organisms with lower map scores probably are more related. Similarly, later you may have multiple organisms and you want to prove that they are all cumulatively related. The computational task is to compute an overall mapping score for a group of strings.

#### Problem Statement: There are K strings Xi from the vocabulary V. Each string Xi has length Ni. Your goal is to map the strings to each other. An easy way to do this is to think of this in two steps – conversion and matching. Conversion is a function F that takes in a string and returns another string. All F(Xi)s have the same length N. N is greater than equal to all Nis. The function F is only allowed to make one change to the original string – it can introduce dashes. It can introduce any number of dashes at any position. The conversion cost of X to F(X) is CC*number of dashes, CC being a constant. Once all strings have been converted the matching step just matches characters at each position. The matching cost between two characters is given by a symmetric function MC(c1, c2) where c1 and c2 are two characters ϵ V U {-}. Matching cost of two strings is the sum of matching costs of their conversions at each position. Finally, the matching cost of K strings is the sum of pairwise matching costs between each pair.

For more detailed problem statement, please go through this [pdf](./A1.pdf)

Note : This program was written to complete the coursework requirements of Artificial Intelligence course [COL333](https://www.cse.iitd.ac.in/~mausam/courses/col333/autumn2023/)
