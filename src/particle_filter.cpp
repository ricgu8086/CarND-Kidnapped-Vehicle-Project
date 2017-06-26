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


/* Prototypes (should be in a .h but the autograder doesn't allow it) */
/* 
 *
 */
inline double squared(double x);

/*****************************************************************************
 * normpdf(X,mu,sigma) computes the probability function at values x using the
 * normal distribution with mean mu and standard deviation std. x, mue and 
 * sigma must be scalar! The parameter std must be positive. 
 * The normal pdf is y=f(x;mu,std)= 1/(std*sqrt(2pi)) e[ -(x−mu)^2 / 2*std^2 ]
 *****************************************************************************/
inline double normpdf(double x, double mu, double std);

inline double normpdfbi(double x, double mu_x, double std_x, \
									double y, double mu_y, double std_y);


/*
 *
 */
vector<double> normalize_vector(vector<double> inputVector);


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

	struct Particle particle_aux;

	for(int i=0; i<num_particles; i++)
	{
		particle_aux.x = dist_x(generator);
		particle_aux.y = dist_y(generator);
		particle_aux.theta = dist_theta(generator);
		particle_aux.weight = 1;
		particle_aux.id = i;

		particles.push_back(particle_aux); // It's added by copy, not by reference. It is OK
	}
	
	weights.resize(num_particles, 1);

	is_initialized = true;	
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
			new_p.x = p.x + velocity * delta_t * cos(p.theta);
			new_p.y = p.y + velocity * delta_t * sin(p.theta);
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

	double dist_min, curr_dist;

	for (auto&& obs: observations)
	{
		dist_min = numeric_limits<double>::max();

		for (auto land: predicted)
		{
			curr_dist = dist(obs.x, obs.y, land.x, land.y);
			
			if(curr_dist < dist_min)
			{
				dist_min = curr_dist;
				obs.id = land.id;
			}
		}
	}
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

	double x, y, mu_x, mu_y, std_x, std_y;
	struct LandmarkObs obs_aux;

	weights.clear();

	// Update: the description of std_landmark is wrong. Reported issue https://github.com/udacity/CarND-Kidnapped-Vehicle-Project/issues/8
	// CORRECT: std_landmark: [x [m], y [m]]
	// WRONG: std_landmark: [standard deviation of range [m], standard deviation of bearing [rad]]
	std_x = std_landmark[0];
	std_y = std_landmark[1];

	for (auto&& p: particles)
	{
		// Transform all observations from vehicle's coordinates to map's coordinates
		vector<LandmarkObs> transformed_obs;

		for (auto&& obs: observations)
		{
			// To do the transformation, first rotate, second translate 
			obs_aux.x = obs.x * cos(p.theta) - obs.y * sin(p.theta) + p.x;
			obs_aux.y = obs.x * sin(p.theta) + obs.y * cos(p.theta) + p.y;
			obs_aux.id = obs.id;

			transformed_obs.push_back(obs_aux);
		}

		// Get all the landmarks within sensor_range distance from the particle (in map's coordinates)
		vector<LandmarkObs> map_land_inrange;

		for(auto map_land: map_landmarks.landmark_list)
		{
			if (dist(p.x, p.y, map_land.x_f, map_land.y_f) <= sensor_range)
			{
				obs_aux.x = map_land.x_f;
				obs_aux.y = map_land.y_f;
				obs_aux.id = map_land.id_i;

				map_land_inrange.push_back(obs_aux);
			}
		}

		// Perform data association with them
		dataAssociation(map_land_inrange, transformed_obs);

		// Update particle weight using a multivariate gaussian taking into account all associated landmarks
		double final_weight = 1;

		vector<int> associations;
		vector<double> sense_x, sense_y;

		for (auto tobs: transformed_obs)
		{
			// Find associated landmark from ID
			for (auto land: map_land_inrange)
				if (tobs.id == land.id)
				{
					mu_x = land.x; // The mean of the Multivariate-Gaussian is the measurement's associated landmark positio
					mu_y = land.y;

					x = tobs.x; // The Multivariate-Gaussian is evaluated at the point of the transformed measurement's position
					y = tobs.y;

					final_weight *= normpdfbi(x, mu_x, std_x, y, mu_y, std_y);

					// For extra information within the simulator
					associations.push_back(tobs.id);
					sense_x.push_back(x);
					sense_y.push_back(y);
				}
		}

		p.weight = final_weight;
		weights.push_back(final_weight);

		SetAssociations(p, associations, sense_x, sense_y);
	}

	weights = normalize_vector(weights);	
}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	mt19937 generator;
    // std::default_random_engine generator; // The previous one was recommended in Udacity code, but I added this one here as a reference (easier to remember)
    discrete_distribution<> distribution(weights.begin(), weights.end());

	vector<Particle> particles2;

	int index;

	for(int i=0; i<weights.size(); i++)
	{
		index = distribution(generator);
		particles2.push_back(particles[index]);
	}

	particles = particles2;
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

vector<double> normalize_vector(std::vector<double> inputVector)
{
	//declare sum:
	double sum = 0.0f;

	//declare and resize output vector:
	vector<double> outputVector;
	outputVector.resize(inputVector.size());

	//estimate the sum:
	for (unsigned int i = 0; i < inputVector.size(); ++i) 
	{
		sum += inputVector[i];
	}

	//normalize with sum:
	for (unsigned int i = 0; i < inputVector.size(); ++i) 
	{
		outputVector[i] = inputVector[i]/sum ;
	}

	//return normalized vector:
	return outputVector ;
}

inline double squared(double x)
{
	return x*x;
}

inline double normpdf(double x, double mu, double std) 
{
	double ONE_OVER_SQRT_2PI = 1/sqrt(2*M_PI) ;

	return (ONE_OVER_SQRT_2PI/std)*exp(-0.5*squared((x-mu)/std));
}

inline double normpdfbi(double x, double mu_x, double std_x, \
										double y, double mu_y, double std_y) 
{
	double ONE_OVER_SQRT_2PI = 1/sqrt(2*M_PI) ;

	return (ONE_OVER_SQRT_2PI/(std_x*std_y))*exp(-0.5*\
			(squared((x-mu_x)/std_x) + squared((y-mu_y)/std_y)));
}