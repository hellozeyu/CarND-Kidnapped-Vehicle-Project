/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <math.h> 

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Number of particles to draw
	num_particles = 100;
	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_thea(theta, std[2]);

	for(int i=0; i<num_particles; ++i) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_thea(gen);
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(1.0);
	}

	// Flag, if filter is initialized
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// noise generation
	default_random_engine gen;
	std::normal_distribution<double> noise_x(0, std_pos[0]);
	std::normal_distribution<double> noise_y(0, std_pos[1]);
	std::normal_distribution<double> noise_theta(0, std_pos[2]);

	for(int i=0; i<num_particles; ++i)
	{
		Particle &p = particles[i];

		double new_x = p.x;
		double new_y = p.y;
		double new_theta = p.theta;

		if(fabs(yaw_rate) < 0.00001){
			new_x = p.x + velocity * delta_t * cos(p.theta);
			new_y = p.y + velocity * delta_t * sin(p.theta);
			new_theta = p.theta;
		} else {
			new_x = p.x + velocity/yaw_rate * (sin(p.theta + (yaw_rate * delta_t)) - sin(p.theta));
			new_y = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + (yaw_rate * delta_t)));
			new_theta = p.theta + yaw_rate * delta_t;
		}

		p.x = new_x + noise_x(gen);
		p.y = new_y + noise_y(gen);
		p.theta  = new_theta + noise_theta(gen);
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	int num_landmarks = map_landmarks.landmark_list.size();
	int num_landmarkObs = observations.size();
	// for every particle
	for (int i = 0; i < num_particles; i++) {
		Particle &particle = particles[i];

		double weight = 1.0;

//        std::vector<LandmarkObs> observations_pred;
		for (int j = 0; j < num_landmarkObs; j++) {
			LandmarkObs obs = observations[j];
			// convert the position to map coordinates
			double map_x = obs.x * cos(particle.theta) - obs.y * sin(particle.theta) + particle.x;
			double map_y = obs.x * sin(particle.theta) + obs.y * cos(particle.theta) + particle.y;

			Map::single_landmark_s closest_landmark;
			double closest_dist;
			// Find the closest landmark to the observation
			for (int map_l = 0; map_l < num_landmarks; map_l++) {
				auto landmark = map_landmarks.landmark_list[map_l];
                // We only want landmarks within the sensor range [x - sensor_range; x+sensor_range] and [y-sensor_range; y+sensor_range]
                // to save some calculation time. Need to tweak the sensor_range parameter in the main.cpp
                if (abs(map_x - landmark.x_f) < sensor_range & abs(map_y - landmark.y_f) < sensor_range){
                    double landmark_dist = dist(map_x, map_y, landmark.x_f, landmark.y_f);
                    if (map_l == 0 || landmark_dist < closest_dist) {
                        closest_landmark = landmark;
                        closest_dist = landmark_dist;
                    }
                }else{
                    continue;
                }

			}
			// Multiply the probability of seeing the closest landmark, with the
			// existing weight, and update the weights
			weight *= (1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1])) *
					  exp(-(
							  (pow(map_x - closest_landmark.x_f, 2.0) / (2.0 * pow(std_landmark[0], 2.0))) +
							  (pow(map_y - closest_landmark.y_f, 2.0) / (2.0 * pow(std_landmark[1], 2.0)))
					  ));
		}
		particle.weight = weight;

	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	std::vector<Particle> resampled_particles;
	discrete_distribution<int> dist_index(weights.begin(), weights.end());

	// Calculate the max weight, and setup the distributions
	double max_weight = 0;
	for (int i=0; i<num_particles; i++) {
		max_weight = max(particles[i].weight, max_weight);
	}
	uniform_real_distribution<double> dist_beta(0, 2.0 * max_weight);

	// Resmaple from lesson 13
	int index = dist_index(gen);
	double beta = 0;
	for (int i=0; i<num_particles; i++) {
		beta += dist_beta(gen);
		while (particles[index].weight < beta) {
			beta = beta - particles[index].weight;
			index = (index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
