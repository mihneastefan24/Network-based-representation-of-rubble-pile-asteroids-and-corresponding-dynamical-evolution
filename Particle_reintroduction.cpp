// ----------------------------------------------------------------------------------
// Code that reintroduces the particles into the GRAINS and creates the desired setup
// ----------------------------------------------------------------------------------
// Martin M.S.
// 2025
// Politecnico di Milano
// -----------------------------------------------------------------------------------
// Chrono Project - https://projectchrono.org/
// -----------------------------------------------------------------------------------

// C++ libraries
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
// GRAINS libraries
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono_irrlicht/ChApiIrr.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/assets/ChColor.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/collision/ChCollisionShapeConvexHull.h"
#include "chrono/collision/ChCollisionShapeSphere.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/collision/ChCollisionSystem.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChVector3.h"
#include "chrono/physics/ChContactMaterialNSC.h"
#include "chrono/physics/ChBody.h"
#include "chrono_vsg/ChVisualSystemVSG.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"

#ifdef  _WIN32
#include <direct.h>
#define chdir _chdir
#else
#include <unistd.h>
#endif //  _WIN32


using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;
using namespace chrono::utils;
using namespace chrono::vsg3d;
using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;


std::vector<double> values_vector(std::string loading_folder) {
    // --------------------------------------------------------------------
    // This function transforms the data from the floder into a C++ vector
    // Input: loading_folder - .txt file path of the desired data
    // Output: data_variable - vector of doubles containing the data from loading_folder
    // --------------------------------------------------------------------

	std::string line;
	std::ifstream data(loading_folder);

	std::getline(data, line);
	std::vector<double> data_variable;

    if (data.is_open()) {
		// Use a while loop together with the getline() function to read the file line by line
		while (std::getline(data, line)) {

			std::stringstream ss(line);
			double value;
			ss >> value;
			if (!ss.fail()) {
				data_variable.push_back(value);
			}
			else {
				std::cerr << "Warning: Invalid data in line: " << line << std::endl;
			}
		}
		// Close the file
		data.close();
		return data_variable;
	}
	else {
		std::cout << "cannot open the file /n";
	}
}

std::vector<std::vector<double>> values_matrix(std::string loading_folder) {
    // --------------------------------------------------------------------
    // This function transforms the data from the floder into a C++ matrix
    // Input: loading_folder - .txt file path of the desired data
    // Output: data_variable - matrix of doubles containing the data from loading_folder
    // --------------------------------------------------------------------

	std::string line;
	std::ifstream data(loading_folder);

	std::getline(data, line);

	std::vector<std::vector<double>> data_variable;

	if (data.is_open()) {
		while (std::getline(data, line))
		{
			std::stringstream ss_pos(line);
			std::vector<double> row;
			double value_pos;

			while (ss_pos >> value_pos)
			{
				row.push_back(value_pos);
			}
			if (!row.empty())
				data_variable.push_back(row);


		}
		data.close();
	}
	else {
		std::cout << "Cannot open the file";
	}
	return data_variable;
}

//     A callback executed at each particle creation can be attached to the emitter.
//     For example, we need that new particles will be bound to Irrlicht visualization:
class MyCreatorForAll : public ChRandomShapeCreator::AddBodyCallback {
public:
	virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
		ChCoordsys<> mcoords,
		ChRandomShapeCreator& mcreator) override {
		// Bind visual model to the visual systemph
		mbody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/stone.png"));
		vis->BindItem(mbody);


		// Bind the collision model to the collision systemph
		if (mbody->GetCollisionModel())
			coll->Add(mbody->GetCollisionModel());

		// Disable gyroscopic forces for increased integrator stability
		mbody->SetUseGyroTorque(false);
	}
	ChVisualSystem* vis;
	ChCollisionSystem* coll;
};


// Define a custom callback to detect contacts
// Custom class that detects the contacts and the contact forces,
// It can also offer the location of each contact w.r.t. the centre  of the particle
class ContactManager : public ChContactContainer::ReportContactCallback {
public:
    // Initialize the function
	ContactManager(ChSystemMulticore* sys) : m_system(sys) {}

	unsigned int GetNumContacts(std::shared_ptr<ChBody> body) const {
		auto search = m_bcontacts.find(body.get());
		return (search == m_bcontacts.end()) ? 0 : search->second;
	}

	std::unordered_set<ChBody*> GetContacts(std::shared_ptr<ChBody> body)
	{
		auto search = m_contactMap.find(body.get());
		if (search != m_contactMap.end())
			return search->second;
		else
			return{};
	}

	std::unordered_set<double> GetContactForces(std::shared_ptr<ChBody> body) {
		auto search = m_forces_restriction.find(body.get());
		if (search != m_forces_restriction.end())
			return search->second;
		else
			return {};
	}

	void Process() {
		m_contactMap.clear();
		m_bcontacts.clear();
		m_forces_restriction.clear();
		std::shared_ptr<ContactManager> shared_this(this, [](ContactManager*) {});
		m_system->GetContactContainer()->ReportAllContacts(shared_this);
	}

private:
	std::unordered_map<ChBody*, std::unordered_set<ChBody*>> m_contact_map;
    // Data that can be extracted from the report
	virtual bool OnReportContact(const ChVector3d& pA, // position of main particle
		const ChVector3d& pB,                          // position particle of contact
		const ChMatrix33<>& plane_coord,               
		const double& distance,                         // Distance
		const double& eff_radius,                       // Effective Radius
		const ChVector3d& cforce,                       // Contact force
		const ChVector3d& ctorque,                      // Contact Torque
		ChContactable* modA,
		ChContactable* modB) override {
		auto bodyA = dynamic_cast<ChBody*>(modA);
		auto bodyB = dynamic_cast<ChBody*>(modB);


		double force_magnitude = cforce.Length();       // abs(F_contact)

        // Impose an error limit 
		if (force_magnitude < 1e-9 || !bodyA || !bodyB)
			return true;

		//	std::cout << "Time: " << m_system->GetChTime() << " Contact between " << bodyA->GetIdentifier()
		//		<< " and " << bodyB->GetIdentifier() << std::endl;
		//	std::cout << "  Contact Force: " << cforce << std::endl;

        // Check if bodyA is in contact with bodyB in the contact report
		auto searchA = m_bcontacts.find(bodyA);
		if (searchA == m_bcontacts.end())
			m_bcontacts.insert(std::make_pair(bodyA, 1));
		else
			searchA->second++;

		auto searchB = m_bcontacts.find(bodyB);
		if (searchB == m_bcontacts.end())
			m_bcontacts.insert(std::make_pair(bodyB, 1));
		else
			searchB->second++;

        // Add the two bodies ids in the maps of each other
		m_contactMap[bodyA].insert(bodyB);
		m_contactMap[bodyB].insert(bodyA);

        // Add the contact forces magnitudes in the maps of each other
		m_forces_restriction[bodyA].insert(force_magnitude);
		m_forces_restriction[bodyB].insert(force_magnitude);

		return true;
	}

	ChSystemMulticore* m_system;
    // Create maps for contact. forces and pairs -> for each particle save the id s of the particle
    // with which it is in contact, and the force magnitude
	std::unordered_map<ChBody*, unsigned int> m_bcontacts;
	std::unordered_map<ChBody*, std::unordered_set<ChBody*>> m_contactMap;
	std::unordered_map<ChBody*, std::unordered_set<double>> m_forces_restriction;
};

ChVector3<> RotateAroundAxis(const ChVector3<>& point, double angle) {
    // ----------------------------------------------------------------------
    // This functions rotate around the z-axis
    // Input: angle: angle of rotation
    //        point: the point of rotation
    // Output: num: the values that is assigned to the parameter in the document  
    // ----------------------------------------------------------------------

	double cosA = cos(angle);
	double sinA = sin(angle);

	// Rotation around Z-axis
	return ChVector3<>(
		cosA * point.x() - sinA * point.y(),  // New X
		sinA * point.x() + cosA * point.y(),  // New Y
		point.z()                             // Z remains unchanged
	);
}

void VideoFrameSave(irr::IrrlichtDevice* device, unsigned int videoframe_each)
{
	unsigned int videoframe_num = 0;
	irr::video::IVideoDriver* driver = device->getVideoDriver();  // Correct way
	if (videoframe_num % videoframe_each == 0) {
		std::filesystem::create_directory("video_capture");
		irr::video::IImage* image = driver->createScreenShot();
		char filename[100];
		sprintf(filename, "video_capture/screenshot%05d.bmp", (videoframe_num + 1) / videoframe_each);
		if (image)
			device->getVideoDriver()->writeImageToFile(image, filename);
		image->drop();
	}
	videoframe_num++;
}


double ReadFFile(std::string parameter, std::string filepath) {
    // ----------------------------------------------------------------------
    // This functions reads the data from a file and assigned the desired value to a double variable
    // Input: parameter: the desired parameter to be looked into the document
    //        filepath:  the path of the document in which to look for parameter
    // Output: num: the values that is assigned to the parameter in the document  
    // ----------------------------------------------------------------------

	std::ifstream thefile(filepath);
	std::string line;
	double num = NULL;
	if (thefile.is_open()) {
		while (getline(thefile, line))
		{
			int l = parameter.length();

			if (parameter == line.substr(0, l))
			{
				std::string numer = line.substr(l + 3, 1000);
				num = stod(numer);
				break;
			}
		}
	}
	else {
		std::cout << filepath << std::endl;
		std::cout << "cannot open file in readfile" << std::endl;
	}
	return num;
}



int main(int argc, char* argv[]) {


	// Initial Data 
	float radius_factor = 10.0;
	float friction = 0.6, rolling_friction = 0.5, spinning_friction = 0.6;
	float cohesion = 0.2, damping = 0.1;
	float complianceT = 0.5, complianceS = 0.5, complianceR = 0.5;
	// Create the Chrono project
	// Declaration of variables

	double restitution = 0.01;

	ChVector3<> center(0, 0, 0);
    // Create the data path through a series for each case, and to read each text file,
    // Not adviceable to put many cases, as the simulation may take a very long time, and fail after
    // several pauses
    //                                                    Name of cases(asteroids)-the name of the folder where each asteroid data has been saved
	std::vector<std::string> asteroid_names_obj = { "Geographos.obj/"};// ,Kleopatra.obj "Arrokoth Stern 2019.obj / ", "Bennu_v20_200k.obj / ", "Geographos Radar - based, mid - res.obj / ", "Apophis Model 1.obj / "};
	// Similarly for the densities, particle shape and the rotation period
    std::vector<std::string> case_names = {"spherical_no_core_monodisperse/"};
	std::vector<std::string> densities_folder = { "rho 3000/", "rho 2400/", "rho 2000/", "rho 1600/", "rho 1200/"}; // , "rho 3000/", "rho 2400/", "rho 2800/"}; "Rotation Period 0.288/","Rotation Period 1.260/", "Rotation Period 1.980/","Rotation Period 3.060/", 
    std::vector<std::string> angvel_folder = { "Rotation Period 3.060/", "Rotation Period 3.240/" }; //  "Rotation Period 0.720/" , "Rotation Period 0.720/", "Rotation Period 1.260/", "Rotation Period 1.980/",
    // Create a folder where to save the name of the data files - Mass, Position, Radius, Orientation etc.
	std::vector<std::string> folder_name = { "Mass.txt","Positions.txt", "Radius.txt" };
	for (auto asteroid_name_obj : asteroid_names_obj)
	{
		for (auto case_name : case_names) {
			std::vector<std::string> asteroid_name = { asteroid_name_obj };
			std::string path_name = "C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/" + case_name + "results/";
            // Create variables in which to save the data
            std::vector<std::vector<double>> positions_apo;
			std::vector<double> mass_apo;
			std::vector<double> radius_apo;

			for (auto dens_folder : densities_folder) {
				for (auto angv_folder : angvel_folder) {
                    // Read the data corresponding to each shape
					std::string path_load = path_name + asteroid_name_obj + dens_folder + angv_folder + "simInputs.txt";
						std::cout << path_load << std::endl;
                    // Extract data from the files
					double G = ReadFFile("Universal gravity constant G", path_load);
					double sf = ReadFFile("Scaling factor ", path_load);
					double sfff = ReadFFile("Scaling factor ^ 3", path_load);
					double density = ReadFFile("Material density", path_load);
					double omega = ReadFFile("Spin rate", path_load);

					for (auto as_name : asteroid_name) {
						for (auto fol_name : folder_name) {
							std::string path = path_name + as_name + fol_name;

                            // Extract the particles'positions
							if (fol_name == "Positions.txt") {
								positions_apo = values_matrix(path);
							}
                            // Extract the particles' masses - unused now
							else if (fol_name == "Mass.txt") {
								mass_apo = values_vector(path);
							}
                            // Extract the particles' radiuses
							else {
								radius_apo = values_vector(path);
							}

						}

                        // Create directories where to save the results and the photos
						std::string  save_directory = "ResultsBackground/" + as_name + angv_folder + dens_folder;
						std::filesystem::create_directories(save_directory);
						std::string save_directory_images = "PhotoBackground/" + as_name + angv_folder + dens_folder;
						std::filesystem::create_directories(save_directory_images);

                        // Here add the data coming from the matlab and that is needed for the simulation the one that adapts the physical characteristics
                        // into the chrono variable
						double OmegaX = 0.0, OmegaY = 0.0, OmegaZ = omega;

                        // Verification
						std::cout << "Number of particles: " << std::endl;
						std::cout << mass_apo.size() << std::endl;
						std::cout << "Total time in seconds: " << std::endl;

                        // Create and set up the multicore system to handle the contacts
						ChSystemMulticoreNSC systemph;
						systemph.SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
						systemph.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
						systemph.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;
						systemph.GetSettings()->solver.max_iteration_normal = 100;
						systemph.GetSettings()->solver.max_iteration_sliding = 100;
						systemph.GetSettings()->solver.max_iteration_spinning = 0;
						systemph.GetSettings()->solver.max_iteration_bilateral = 150;
						systemph.GetSettings()->solver.tolerance = 1e-3;

						systemph.SetNumThreads(6);
						systemph.SetTimestepperType(chrono::ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
						systemph.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

                        // Create the contact material
						auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
						particle_mat->SetFriction(0.6f);
						particle_mat->SetRollingFriction(rolling_friction);
						particle_mat->SetSpinningFriction(spinning_friction);
						particle_mat->SetComplianceRolling(complianceR);
						particle_mat->SetComplianceSpinning(complianceS);

						auto  OmegaSystem = ChVector3d(OmegaX, OmegaY, OmegaZ);

                        // Create the particles and attach them to the system
						for (size_t i = 0; i < positions_apo.size(); i++) {
							auto sphere = chrono_types::make_shared<chrono::ChBodyEasySphere>(radius_apo[i], density, true, true, particle_mat);
							ChVector3<> positions = { positions_apo[i][0], positions_apo[i][1], positions_apo[i][2] };
							sphere->SetPos(positions);
							sphere->SetPosDt({ 0.0,0.0,0.0 });
							sphere->SetPosDt2({ 0,0,0 });
							sphere->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/stone.png"));
							systemph.AddBody(sphere);
						}

						auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
						vis->AttachSystem(&systemph);
						vis->SetWindowSize(800, 600);
						vis->SetWindowTitle("NSC Multicore No Core");
						vis->Initialize();
						vis->AddCamera(ChVector3d(50,50, 50));
						vis->AddTypicalLights();
						vis->AddSkyBox(GetChronoDataFile("textures/white.jpg"));
		

                        // Initialize the contact manager - it is important in order to extract
                        // the contact pairs and contact forces
						ContactManager manager(&systemph);

                        // Create the documents in which to save the data
	.   
						std::string path = save_directory + "Positions.txt"; std::ofstream positions(path);
						path = save_directory + "Velocities.txt";		     std::ofstream velocity(path);
						path = save_directory + "AngularMomentum.txt";	     std::ofstream angular(path);
						path = save_directory + "ContactForces.txt";	     std::ofstream contactf(path);
						path = save_directory + "ContactPositions.txt";	     std::ofstream contactp(path);
						path = save_directory + "AccumulatedForce.txt";	     std::ofstream accumulatedforce(path);
						path = save_directory + "NrContacts.txt";	         std::ofstream nrcontacts(path);
						path = save_directory + "Radius.txt";	             std::ofstream radius(path);
						path = save_directory + "Mass.txt";	                 std::ofstream mass(path);
						path = save_directory + "Inertia.txt";			     std::ofstream inertia(path);
						path = save_directory + "MyContactsNr.txt";		     std::ofstream mycontacts(path);
						path = save_directory + "EachForce.txt";		     std::ofstream eachforce(path);
						path = save_directory + "AngularVelocity.txt";	     std::ofstream angularvel(path);
						path = save_directory + "ContactPairs.txt";	         std::ofstream pairs(path);

                        //Simulation parameter
                        int timestep = 0.1;
						int time = 0;
						int max_time = 3600 * 8; // Set maximum time of the simulation
						std::cout << max_time << std::endl;
						double G_constant = G;

						// Simulation loop
						double screenshot_interval = 3600.0; // Every 3600 simulation seconds
						double next_screenshot_time = 0.0;
						int out_time = 0;
						double out_step = 1;
						unsigned int iterafter = 0;
						double max_dist = 0;

						int frame = 0;
						double contact = 0;
						int timet = 0;
						bool impl = true;
						ChVector3 camera_position;

                        // Impose a limit of 8 hours for each simulation -> max_time
						while (vis->Run() && systemph.GetChTime() <= max_time) {

							vis->BeginScene();
							vis->Render();
							vis->EndScene();

							// Apply custom forcefield (brute force approach..)
							// A) reset 'user forces accumulators':
							unsigned int iter = 0;
							for (auto body : systemph.GetBodies()) {
								body->EmptyAccumulators();
								//mass << iter << "\t" << body->GetMass() << std::endl;
								iter++;
							}

							// B) store user computed force: - brute force approach
							auto& bodies = systemph.GetBodies();
							size_t N = bodies.size();
							std::vector<ChVector3d> force_buffers(N, VNULL);
							for (size_t i = 0; i < N; i++) {
								auto abodyA = bodies[i];
								for (size_t j = i + 1; j < N; j++) {
									auto abodyB = bodies[j];
									ChVector3d D_attract = abodyB->GetPos() - abodyA->GetPos();
									double r_attract = 1 / D_attract.Length();
									double f_attract = G_constant * (abodyA->GetMass() * abodyB->GetMass()) * r_attract * r_attract;

									ChVector3d F_attract = (D_attract * r_attract) * f_attract;

									abodyA->AccumulateForce(F_attract, abodyA->GetPos(), false);
									abodyB->AccumulateForce(-F_attract, abodyB->GetPos(), false);


                                    // Calculate the maximum distance in the system to have a continuosly
                                    // updating camera position in the visualization system
									max_dist = (max_dist < std::sqrt(std::pow(r_attract, 2))) ? std::sqrt(std::pow(r_attract, 2)) : max_dist;
									if (max_dist < 150) {
										max_dist = max_dist;
									}
									else if (max_dist < 200) {
										max_dist = 150;
									}
									else if (max_dist < 250) {
										max_dist = 200;
									}
									else if (max_dist < 300) {
										max_dist = 250;
									}
									else if (max_dist < 350) {
										max_dist = 300;
									}
									else
									{
										max_dist = 350;
									}
                                    // Setup maximum limiting to be able to see the particles
									if (max_dist > 350)
										max_dist = 350;
								}
							}

							// Perform the integration timestep
							systemph.DoStepDynamics(timestep);

							double time = systemph.GetChTime();

							// Get the forces and the contact pairs
							manager.Process();

                            // Extract the forces and the contact pairs every 300 seconds
							if (round(100 * fmod(systemph.GetChTime(), 300)) / 100.0 == 0.0) {

                                for (auto body : systemph.GetBodies())
								{
									auto bodyContacts = manager.GetContacts(body);
									auto forceContacts = manager.GetContactForces(body);

									pairs << systemph.GetChTime() << "\t" << body->GetIndex() << "\t";
									eachforce << systemph.GetChTime() << "\t" << body->GetIndex() << "\t";
									for (auto contact : bodyContacts) {
										pairs << contact->GetIndex() << "\t";
									}

									for (const auto& force : forceContacts)
										eachforce << force << "\t";
									eachforce << "\n";
									pairs << "\n";

								}
							}

                            // Set the asteroid rotation after 50 secs to allow the particles to come all in contact
							if (systemph.GetChTime() > 50 && impl) {
								impl = false;
								for (auto& body : systemph.GetBodies()) {
									ChVector3d r = body->GetPos();
									ChVector3d v_init = Vcross(ChVector3d(OmegaSystem), r);
									body->SetPosDt(v_init);
								}

							}
                            // Save the photos every hour or whenever desired through the next_screenshot_time
							if (time >= next_screenshot_time) {
								vis->WriteImageToFile(save_directory_images + "simulation_frame_" + std::to_string((int)(time)) + ".png");
								next_screenshot_time += screenshot_interval;  // schedule next screenshot
								camera_position = vis->GetCameraPosition();
								camera_position += 20;
								vis->SetCameraPosition(ChVector3d(max_dist - 5, max_dist - 5, max_dist - 5));
							}

                            // Save the desired data every 60 time steps 
							if (round(100 * fmod(systemph.GetChTime(), 300)) / 100.0 == 0.0)
							{
								systemph.CalculateContactForces();
								for (auto body : systemph.GetBodies()) {
									real3 frc = systemph.GetBodyContactForce(body);

									positions << body->GetPos().x() << "\t" << body->GetPos().y() << "\t" << body->GetPos().z() << "\t";
									velocity << body->GetPosDt().x() << "\t" << body->GetPosDt().y() << "\t" << body->GetPosDt().z() << "\t";
									angular << body->GetAngVelLocal().x() << "\t" << body->GetAngVelLocal().y() << "\t" << body->GetAngVelLocal().z() << "\t";
								//	contactf << body->GetContactForce().x() << "\t" << body->GetContactForce().y() << "\t" << body->GetContactForce().z() << "\t";
									contactf << frc.x << "\t" << frc.y << "\t" << frc.z << "\t";
									accumulatedforce << body->GetAccumulatedForce().x() << "\t" << body->GetAccumulatedForce().y() << "\t" << body->GetAccumulatedForce().z() << "\t";
									mass << body->GetMass() << "\t";
									inertia << body->GetInertiaXX().x() << "\t" << body->GetInertiaXX().y() << "\t" << body->GetInertiaXX().z() << "\t";
									angularvel << body->GetAngVelParent().x() << "\t" << body->GetAngVelParent().y() << "\t" << body->GetAngVelParent().z() << "\t";
								}

								for (int i = 0; i < systemph.GetBodies().size(); i++)
								{		radius << radius_apo[i] << "\t";
								}
								
								mycontacts << contact << "\n";
								nrcontacts << systemph.GetNumContacts() << "\n";
								mass << "\n";
								radius << "\n";
								positions << "\n";
								velocity << "\n";
								angular << "\n";
								contactf << "\n";
								accumulatedforce << "\n";
								inertia << "\n";
								angularvel << "\n";
								std::cout << systemph.GetChTime() << "\t" << time << std::endl;
							}
                            // Go to next time step
							timet += timestep;
							frame++;

							if (time >= max_time)  // Check if the visualization is still running
							{
								std::cout << "MAX TIME" << "\t" << path_name + as_name << std::endl;
								vis->Quit();
							}
							
						}
						positions.close();
						velocity.close();
						angular.close();
						contactf.close();
						contactp.close();
						accumulatedforce.close();
						nrcontacts.close();
						radius.close();
						mass.close();
						inertia.close();
						mycontacts.close();
						eachforce.close();
						angularvel.close();
                        // Important to move in the next simulation of the for, all bodies must be cleared
                        // and the system cleared
						systemph.Clear();
						systemph.RemoveAllBodies();
					}
				}
			}
		}
	}
	return 0;
}