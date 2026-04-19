// -------------------------------------------------------------------------------
// Code that generates particle and let them aggregate into a rubble-pile
// -------------------------------------------------------------------------------
// Martin M.S.
// 2025
// Politecnico di Milano
// -------------------------------------------------------------------------------
// Chrono Project - https://projectchrono.org/
// -------------------------------------------------------------------------------

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono_vsg/ChVisualSystemVSG.h"


#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>
#include "set"
#include "string"
#include "random"
#include "vector"

using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;
using namespace chrono::utils;
using namespace chrono::vsg3d;

// Define the global variables
double density = 2800;
double restitution = 0.1;
double adhesion = 0.6;
double adhesionMult = 0.0;
double Kn = 100000;
double Kt = 100000;
double Gn = 40;
double Gt = 20;

class ContactBodyReporter : public chrono::ChContactContainer::ReportContactCallback {
public:
    virtual bool OnReportContact(
        const ChVector3<>& pA,            // contact point on object A
        const ChVector3<>& pB,            // contact point on object B
        const ChMatrix33<>& plane_coord, // contact plane coords (normal, U, V)
        const double& distance,          // penetration distance
        const double& eff_radius,        // effective radius of curvature
        const ChVector3<>& react_forces,  // forces in contact plane
        const ChVector3<>& react_torques, // torques in contact plane
        chrono::ChContactable* objA,     // contactable object A
        chrono::ChContactable* objB      // contactable object B
    ) override {
        // Attempt to cast contactable objects to ChBody
        auto bodyA = dynamic_cast<chrono::ChBody*>(objA);
        auto bodyB = dynamic_cast<chrono::ChBody*>(objB);

        if (bodyA && bodyB) {
            std::cout << "Contact between Body A and Body B:" << std::endl;
            std::cout << " - Body A ID: " << bodyA->GetIdentifier() << std::endl;
            std::cout << " - Body B ID: " << bodyB->GetIdentifier() << std::endl;
        }
        else {
            std::cout << "Contact involves non-body objects." << std::endl;
        }
        return true; // Continue reporting contacts
    }
};



void AddParticles(ChSystem* sys,int& Nbodies,  bool cont_case) {

    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
    particle_mat->SetRestitution(restitution);
//    particle_mat->SetAdhesion(adhesion);
//   particle_mat->SetKn(Kn);
//   particle_mat->SetKt(Kt);
//   particle_mat->SetGn(Gn);
//    particle_mat->SetGt(Gt);

    for (int i = 0; i < Nbodies; i++)
    {
        auto body = chrono_types::make_shared<ChBodyEasySphere>(2, density, true, true, particle_mat);
        body->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
        double x = 150 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
        double y = 150 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
        double z = 150 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
        if (fabs(x) < 60 && cont_case) {
            if (x < 0)
                x -= 70;
            else
                x += 70;
        }
        if (y < 70 && cont_case) {
            if (y < 0)
                y -= 70;
            else
                y += 70;
        }
        if (z < 70 && cont_case) {
            if (z < 0)
                z -= 70;
            else
                z += 70;
        }
        body->SetPos(ChVector3d(x, y, z));
        sys->AddBody(body);
    }
}

void AddParticlesCont(ChSystem* sys, std::vector<std::vector<double>> pos, std::vector<std::vector<double>> vel) {

    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
    particle_mat->SetRestitution(restitution);
//    particle_mat->SetAdhesion(adhesion);
//    particle_mat->SetKn(Kn);
//    particle_mat->SetKt(Kt);
//    particle_mat->SetGn(Gn);
//    particle_mat->SetGt(Gt);

    for (int i = 0; i < 2500; i++)
    {
        auto body = chrono_types::make_shared<ChBodyEasySphere>(2, density, true, true, particle_mat);
        body->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
        body->SetPos(ChVector3d(pos[i][1], pos[i][2], pos[i][3]));
        body->SetPosDt(ChVector3d(vel[i][1], vel[i][2], vel[i][3]));
        sys->AddBody(body);
    }
}



int main(int argc, char* argv[]) {

    // Create a Chrono physical system
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetNumThreads(8);
    
  /*  sys.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    sys.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::PRIMS;
    sys.GetSettings()->solver.max_iteration_normal = 50;
    sys.GetSettings()->solver.max_iteration_sliding = 50;
    sys.GetSettings()->solver.max_iteration_spinning = 0;
    sys.GetSettings()->solver.max_iteration_bilateral = 50;
    sys.GetSettings()->solver.tolerance = 1e-3;
*/
    // Set the systems data
    int cont = 0;
    int Nbodies = 10000 ; 

    
    AddParticles(&sys, Nbodies,false);

    // Write the documents
    std::ofstream pos("Results/Position_Fixed_tstep_small_NSC.txt");
    std::ofstream vel("Results/Velocity_Fixed_tstep_small_NSC.txt");
    std::ofstream mass("Results/Mass_Fixed_tstep_small_NSC.txt");

    // Simulation loop setup
    double timestep = 0.5;
    // If photos are wanted
    double screenshot_interval = 3600.0; // Every 3600 simulation seconds
    double next_screenshot_time = 0.0;
    int out_time = 0;
    double out_step = 100;
    double G_constant = 6.6743e-11;  // gravitational constant 
    unsigned int iterafter = 0;

    // Create the VSG visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemVSG>();
    vis->AttachSystem(&sys);
    vis->SetWindowTitle("Aggregate NSC");
    vis->AddCamera(ChVector3d(0, 150, -150));
    vis->SetWindowSize(ChVector2i(800, 600));
    vis->SetWindowPosition(ChVector2i(100, 100));
    vis->SetClearColor(ChColor(0.8f, 0.85f, 0.9f));
    vis->SetUseSkyBox(true);  // use built-in path
    vis->SetCameraVertical(CameraVerticalDir::Y);
    vis->SetCameraAngleDeg(40.0);
    vis->SetLightIntensity(1.0f);
    vis->SetLightDirection(1.5 * CH_PI_2, CH_PI_4);
    vis->SetShadows(true);
    vis->SetWireFrameMode(false);
    vis->Initialize();

    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    std::cout << sys.GetBodies().size() << std::endl;
    sys.SetSolverType(ChSolver::Type::PARDISO_MKL);
    sys.GetSolver()->AsIterative()->SetMaxIterations(100);

    // while loop initialization
    unsigned int Tstep_k = 0;
    int hours = 0;
    while (vis->Run()) {

        vis->BeginScene();
        vis->Render();
        vis->EndScene();


        unsigned int iter = 0;
        for (auto body : sys.GetBodies()) {
            body->EmptyAccumulators();
            iter++;
        }
            // Create gravitational attraction brute-force
        for (unsigned int i = 0; i < sys.GetBodies().size(); i++) {
            auto abodyA = sys.GetBodies()[i];
            for (unsigned int j = i + 1; j < sys.GetBodies().size(); j++) {
                auto abodyB = sys.GetBodies()[j];
                ChVector3d D_attract = abodyB->GetPos() - abodyA->GetPos();
                double r_attract = D_attract.Length();
                double f_attract = G_constant * (abodyA->GetMass() * abodyB->GetMass()) / (std::pow(r_attract, 2));
                ChVector3d F_attract = (D_attract / r_attract) * f_attract;

                abodyA->AccumulateForce(F_attract, abodyA->GetPos(), false);
                abodyB->AccumulateForce(-F_attract, abodyB->GetPos(), false);
            }
        }

        // Perform the integration timestep
        sys.DoStepDynamics(timestep);

        double time = sys.GetChTime();
        // Save the data every hour - this deletes the previous hour data savings-> for possible restarts
        if (int(time) % 3600 == hours)
        {
            std::ofstream pos6("Position_Fixed_tstep_small_NSC_" + std::to_string(int(time / 3600)) + ".txt");
            std::ofstream vel6("Velocity_Fixed_tstep_small_NSC_" + std::to_string(int(time / 3600)) + ".txt");
            for (auto body : sys.GetBodies()) {
                //  volume = density * body->GetMass();
                pos6 << iterafter << "\t" << body->GetPos().x() << "\t" << body->GetPos().y() << "\t" << body->GetPos().z() << std::endl;
                vel6 << iterafter << "\t" << body->GetPosDt().x() << "\t" << body->GetPosDt().y() << "\t" << body->GetPosDt().z() << std::endl;
                iterafter++;
            }
           pos6.close();
           vel6.close();
           hours++;
        }
    }
    
    // Take the values of the particles at the end when the body is aggregated
    iterafter = 0;
    for (auto body : sys.GetBodies()) {
        body->EmptyAccumulators();
        //  volume = density * body->GetMass();
        
        pos << iterafter << "\t" << body->GetPos().x() << "\t" << body->GetPos().y() << "\t" << body->GetPos().z() << std::endl;
        vel << iterafter << "\t" << body->GetPosDt().x() << "\t" << body->GetPosDt().y() << "\t" << body->GetPosDt().z() << std::endl;
        mass << iterafter << "\t" << body->GetMass() << std::endl;
        iterafter++;
        
    }
    mass.close();
    pos.close();
    vel.close();
    return 0;
}
