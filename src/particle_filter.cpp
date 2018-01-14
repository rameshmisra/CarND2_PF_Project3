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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	num_particles = 25;
	
	// Create normal distributions:
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
    //cout<< "x y theta gps values " << x<<" "<< y<< " "<<theta<< endl;

	for (int i = 0; i < num_particles; i++) 
	{ 
		Particle P;

		P.id = i;
		P.x = dist_x(gen);
		P.y = dist_y(gen);
		P.theta = dist_theta(gen);
		P.weight = 1.0;
        //cout<< "particles: " <<P.x << " " << P.y << " " <<P.theta << " " << P.weight<< endl;
		particles.push_back(P);
		weights.push_back(1.0);
	}

	is_initialized = true;
	

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	

	//vector<Particle> particles = pf.particles;
	//num_particles = particles.size();
	double x_new;
	double y_new;
	double theta_new;

	for (int i = 0; i < num_particles; i++)
	{ 
		if (fabs(yaw_rate) > 0.001)
		{
			x_new = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_new = particles[i].y - (velocity/yaw_rate)*(cos(particles[i].theta + yaw_rate*delta_t) - cos(particles[i].theta));
			theta_new = particles[i].theta + yaw_rate*delta_t;
		}
		else
		{
			x_new = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			y_new = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			theta_new = particles[i].theta;
		}

		//add noise
		normal_distribution<double> dist_x(x_new, std_pos[0]);
	    normal_distribution<double> dist_y(y_new, std_pos[1]);
	    normal_distribution<double> dist_theta(theta_new, std_pos[2]);

	    particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		//cout<< "moved particles: "<<particles[i].x <<" "<< particles[i].y <<" "<<particles[i].theta <<" "<< particles[i].weight<< endl;

    }

    //cout<< "pf.prediction completed " << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	
    double min_dist = 50*20.0;   //sensor_range = 50; start with large distance
    double distance;
    int associated_index = 99;

    for (int i=0; i< observations.size(); i++)
    {
    	for (int j=0; j < predicted.size(); j++)
		{
			distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if(distance <= min_dist)
			{
				min_dist = distance;
				associated_index = predicted[j].id;	
			}
        
		}
		observations[i].id = associated_index;
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

    std::vector<LandmarkObs> XformedObs;   //vector of transformed obs in map coords using struct of LandmarkObs
    std::vector<LandmarkObs> withinrangeLMs;
    double denominator = 2.0*M_PI*std_landmark[0] * std_landmark[1];  //for weight calculation

	for (int i=0; i < num_particles; i++)
	{
		double x_particle = particles[i].x;
        double y_particle = particles[i].y;
        double theta_particle = particles[i].theta;

		
		withinrangeLMs.clear();
		

		for (int k=0; k< map_landmarks.landmark_list.size(); k++)
		{
			if (dist(x_particle, y_particle, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f) < sensor_range)
			{
				withinrangeLMs.push_back(LandmarkObs{map_landmarks.landmark_list[k].id_i, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f});
			//cout<<"Px&y "<<x_particle<<" "<< y_particle<<" "<<map_landmarks.landmark_list[k].id_i<<endl;
			}
		}
		//cout<< "Total landmarks "<< map_landmarks.landmark_list.size()<<" LMs in range "<<withinrangeLMs.size()<< endl;

		XformedObs.clear();
		LandmarkObs obs_new;  //vector of landmark observartions for particle i.
		particles[i].weight = 1.0;

		for (int j=0; j < observations.size(); j++)
		{
			double x_obs = observations[j].x;
			double y_obs = observations[j].y;

			obs_new.id = observations[j].id;
			obs_new.x = x_particle + x_obs * cos(theta_particle) - y_obs * sin(theta_particle);
			obs_new.y = y_particle + x_obs * sin(theta_particle) + y_obs * cos(theta_particle);

			XformedObs.push_back(obs_new);
			//cout<< "Obs Landmarks: "<< x_obs<<","<<y_obs<<" xformed: " <<obs_new.x << ","<<obs_new.y<<" "<< x_particle<< endl;


			double min_dist = sensor_range; 
            double distance;
            int associated_index = 99;
            double lx;
            double ly;

			for (int n=0; n < withinrangeLMs.size(); n++)
		    {
			    distance = dist(obs_new.x, obs_new.y, withinrangeLMs[n].x, withinrangeLMs[n].y);
				if(distance <= min_dist)
				{
					min_dist = distance;
					associated_index = withinrangeLMs[n].id;
					lx = withinrangeLMs[n].x;
					ly = withinrangeLMs[n].y;
					//cout<< "associated index: "<< associated_index<< "  min dist: "<< min_dist	<<endl;
				}
		    }
		    //cout<< "Final associated index: "<< associated_index<< "  min dist: "<< min_dist	<<endl;

		    //probability & weight calculation

		    //cout<<"Associated Landmark x,y "<<lx<<", "<<ly<<endl;

            double xterm = (obs_new.x - lx)/std_landmark[0];
            double yterm = (obs_new.y - ly)/std_landmark[1];
            //cout<< "xterm "<< xterm <<"yterm "<<yterm<< endl;
            //cout<< "obs new x,y "<< obs_new.x<<"  "<< obs_new.y <<endl;

		    particles[i].weight *= exp(-0.5*(xterm*xterm + yterm*yterm))/denominator;
		    weights[i] = particles[i].weight;
		}

		//dataAssociation(withinrangeLMs, XformedObs); 
		

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resampled_particles;
	std::random_device rd;
    std::mt19937 gen(rd());      //19937 is a mersenne prime seed # used to generate random number
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    //std::map<int, int> m;
    for (int i=0; i<num_particles; i++) 
    {
        resampled_particles.push_back(particles[dist(gen)]);
    }
    particles = resampled_particles;

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
