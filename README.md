# Network-based-representation-of-rubble-pile-asteroids-and-corresponding-dynamical-evolution
The codes included in this repository are the codes used for analysing the applications of network theory into gravitational aggregates, and instructions for their use.

# Codes Descriptions

Initial_aggregate.cpp - creates an aggregate of 10 000 particles, which are allowed to interact through the brute force until a stable gravitational-aggregate is formed. It has no required input, and creates as an output .txt files that contains the values characterising each particle in the final state of the aggregate. The user must allow the simulation to run until it observes that the particles have merged.

Shape_and_Sizing.m    - loads the data obtained from the previous code and shapes the gravitational aggregate into the shape of the wanted asteroid - in this case Kleopatra, Apophis, Geographos, Bennu, Arrokoth. Once the asteroid has been shaped, the specific data for each asteroid is calculated in order to impose a global total volume to each asteroid, which allows for an accurate comparison between the differen cases in the following steps.
The parameters selected in this code have been selected based on the literature survey, thus they can be modified if other data or other simulations are of interest to the user.

Particle_reintroduction.cpp - takes the data obtained in the previously two steps and creates a new GRAINS systems in which the desired density and angular velocity is imposed. The particles' state is taken based on the previous two codes, and the granular data is considered to be based on literature and is kept constant with the data set in Initial_aggregate.cpp.

Graph_creation_analysis.m - this codes uploads the data obtained in the previous step and creates the graph, analysing the specific network metrics desired. Others methods can be analysed that are already implemented or are easy to implement in matlab.

Edge_elimination_process.m - along with the Graph_creation_analysis.m, it simulates the edge elimination process, and is a continuation of the previous code. 


NOTE: In order to run the C++ codes, one must install the Chrono Project ("projectchrono.org") and install as instructed the libraries for the simulation. Furthermore, other libraries that are based on GRAINS should be considered ("Ferrari, F., Tasora, A., Masarati, P., & Lavagna, M. 2017, Multibody System Dynamics, 39, 3").
