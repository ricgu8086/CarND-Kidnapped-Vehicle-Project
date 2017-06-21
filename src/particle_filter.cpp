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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 50; // TODO parameter to be tuned after some testing

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	
	default_random_engine generator;

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	struct particle aux;

	for(int i=0; i<num_particles; i++)
	{
		particle.x = dist_x(generator);
		particle.y = dist_y(generator);
		particle.theta = dist_theta(generator);
		particle.weight = 1;

		particles.push_back(particle); // It's added by copy, not by reference. It is OK
	}
	
	weights.resize(num_particles, 1);
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	default_random_engine generator;

	struct Particle new_p;

	for(auto&& p: particles)
	{

		if (fabs(yaw_rate) < 1e-6) // We can consider yaw_rate = 0
		{
			new_p.x = p.x + v * delta_t * cos(p.theta);
			new_p.y = = p.y + v * delta_t * sin(p.theta);
			new_p.theta = p.theta;
		}
		else
		{
			new_p.x = p.x + velocity/yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
			new_p.y = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
			new_p.theta = p.theta + yaw_rate * delta_t;
		}

		normal_distribution<double> dist_x(new_p.x, std_x);
		normal_distribution<double> dist_y(new_p.y, std_y);
		normal_distribution<double> dist_theta(new_p.theta, std_theta);

		p.x = dist_x(generator);
		p.y = dist_y(generator);
		p.theta = dist_theta(generator);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) 
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) 
{
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
}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	mt19937 generator;
    // std::default_random_engine generator; // The previous one was recommended in Udacity code, but I added this one here as a reference (easier to remember)
    discrete_distribution<int> distribution(weights);

	vector<double> weights2;

	int index;

	for(int i=0; i<weights.size(); i++)
	{
		index = distribution(generator);
		weights2.push_back(weights[index]);
	}

	weights = weights2;
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
