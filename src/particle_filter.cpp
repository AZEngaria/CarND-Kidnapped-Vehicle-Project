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
#include <float.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
const double pi= 3.14159265358979323846;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	cout<<"Initialization: "<<endl;
	default_random_engine gen;
	
	//The following lines creates a normal (Gaussian) distribution for x, y  and theta.
	normal_distribution<double> dist_x(x, std[0]);	
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	num_particles = 100;
	
	Particle p;
	
	for(int i=0; i < num_particles; i++){		
		p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
		//p.theta = p.theta % (2 * pi);
		p.id = i;
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(1.0);
	}
	is_initialized = true;
	cout<<"Particles Initialized "<<endl;
	return;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	
	cout<<"Predictions: "<<endl;
	
	for(int i=0; i < num_particles; i++){
		Particle p = particles[i];
		double x, y, theta;
		//Motion Upate
		if(fabs(yaw_rate) > 0.0001){  
			x = p.x + (velocity/yaw_rate) * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
			y = p.y + (velocity/yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
			theta = p.theta + yaw_rate * delta_t;
		} else {
			x = p.x + velocity * delta_t * cos(p.theta);
			y = p.y + velocity * delta_t * sin(p.theta);
			theta = p.theta;
		}
		
		
		//The following lines creates a normal (Gaussian) distribution for x, y  and theta.
		normal_distribution<double> dist_x(x, std_pos[0]);	
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);
		
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		//p.theta = p.theta % ( 2*pi);
		
	}
	return;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
		
	cout<<"Data Assosciation: "<<endl;
	int index;
	for(int i=0; i < predicted.size(); i++){
		double min = DBL_MAX;
		for(int j=0; j < observations.size(); j++){
			double d = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
			if(d < min){
				min = d;
				index = j;
			}
		}
		observations[index].id = predicted[i].id;
		
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	// Converting to Map coordinate system
	cout<<"Update Weights: "<<endl;
	
	for(int i=0; i < num_particles; i++){
		Particle p = particles[i];
		double theta = p.theta;
		std::vector<LandmarkObs> Tobservations;
		LandmarkObs transformedLandmark;
		for(int j = 0; j < observations.size() ; j++){
			LandmarkObs observedlandmark = observations[j];
			double x = p.x + observedlandmark.x * cos(theta) - sin(theta) * observedlandmark.y;
			double y = p.y + observedlandmark.x * sin(theta) + cos(theta) * observedlandmark.y;
			transformedLandmark.x = x;
			transformedLandmark.y = y;
			Tobservations.push_back(transformedLandmark);
		}	
		
		std::vector<LandmarkObs> MapLandmarks;		
		LandmarkObs MapLandmark;
		
		for(int j=0; j< map_landmarks.landmark_list.size(); j++){
			MapLandmark.x = map_landmarks.landmark_list[j].x_f;
			MapLandmark.y = map_landmarks.landmark_list[j].y_f;
			MapLandmark.id = map_landmarks.landmark_list[j].id_i;
			if(dist(p.x, p.y, MapLandmark.x, MapLandmark.y) <= sensor_range)
				MapLandmarks.push_back(MapLandmark);
		}
		
		
		dataAssociation(MapLandmarks,Tobservations);
		
		double sig_x= std_landmark[0];
		double sig_y= std_landmark[1];
		float x_obs;
		float y_obs;
		float mu_x;
		float mu_y;
		double final_weight = 1;
		double weight;
		double gauss_norm;
		double exponent;
		for(int k=0; k < Tobservations.size(); k++){
			for(int j = 0; j < MapLandmarks.size(); j++){
				if(Tobservations[k].id == MapLandmarks[j].id){
					x_obs = Tobservations[k].x;
					y_obs = Tobservations[k].y;
					mu_x = MapLandmarks[j].x;
					mu_y = MapLandmarks[j].y;
					cout<<"Found" << endl;
					break;
				}
			}
			//calculate normalization term
			gauss_norm = (1/(2 * pi * sig_x * sig_y));

			//calculate exponent
			exponent = (pow((x_obs - mu_x),2)/(2 * pow(sig_x,2))) + (pow((y_obs - mu_y),2)/(2 * pow(sig_y,2)));

			//calculate weight using normalization terms and exponent
			weight = gauss_norm * exp(-exponent);

			final_weight *= weight;
		}
		particles[i].weight = final_weight;
		weights.push_back(final_weight);
	}
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	cout<<"Resample: "<<endl;
	
	std::default_random_engine generator;
	std::discrete_distribution<> distribution (weights.begin(), weights.end());	
	
	std::vector<Particle> new_particles;
	
	for (int i=0; i<num_particles; ++i) {		
		int number = distribution(generator);
		Particle p = particles[number];
		new_particles.push_back(p);
	}
	particles = move(new_particles);

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
