# Network-based-representation-of-rubble-pile-asteroids-and-corresponding-dynamical-evolution
The codes included in this repository are the codes used for analysing the applications of network theory into network theory, and instructions for their use.

# Codes Descriptions

Initial_aggregate.cpp - creates an aggregate of 10 000 particles, which are allowed to interact through the brute force until a stable gravitational-aggregate is formed.

Shape_and_Sizing.m    - loads the data obtained from the previous code and shapes the gravitational aggregate into the shape of the wanted asteroid - in this case Kleopatra, Apophis, Geographos, Bennu, Arrokoth. Once the asteroid has been shaped, the specific data for each asteroid is calculated in order to impose a global total volume to each asteroid, which allows for an accurate comparison between the differen cases in the following steps.

Particle_reintroduction.cpp - takes the data obtained in the previously two steps and creates a new GRAINS systems in which the desired density and angular velocity is imposed.

Graph_creation_analysis.m - this codes uploads the data obtained in the previous step and creates the graph, analysing the specific network metrics desired

Edge_elimination_process.m - along with the Graph_creation_analysis.m, it simulates the edge elimination process.
