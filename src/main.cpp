#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

/* *********************************************************************
 * 
 *                     Trajectory Helper Functions 
 * 
 * ********************************************************************/
 
 // For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Elucidian distance
double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Get the nearest waypoint ignoring if it's behind or in front of the vehicle
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	
	return closestWaypoint;
	
}

// Get the nearest waypoint in front of the vehicle
int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Get two reference points that can be used for drawing the trajectory near to the car's current position
vector<vector<double>> get_starter_reference_points(const vector<double>& prev_path_x, const vector<double>& prev_path_y, const double& car_x, const double& car_y, const double& car_yaw) {
	
	// Create placeholder vectors for storing reference points to use for guiding the trajectory
	vector<double> ref_points_x;
	vector<double> ref_points_y;
	
	// Check to see if we have previous path points left over... i.e. not the first cycle
	if (min(prev_path_x.size(), prev_path_y.size()) >= 2) {
		
		// If so, use the soonest two trajectory points to define the new trajectory's start
		ref_points_x = prev_path_x;
		ref_points_y = prev_path_y;
		
	} else {
		
		// If we don't have a previous path to work with, just use the current position and yaw
		ref_points_x.push_back(car_x - cos(car_yaw));
		ref_points_x.push_back(car_x);
		ref_points_y.push_back(car_y - sin(car_yaw));
		ref_points_y.push_back(car_y);
		
	}
	
	return {ref_points_x, ref_points_y};
	
}

// Use spline to create an xy path for a given lane
tk::spline get_spline(vector<vector<double>> ref_points, int ref_lane, const double& car_x, const double& car_y, const double& car_s, const double& car_yaw, const std::vector<double>& map_x, const std::vector<double>& map_y, const std::vector<double>& map_s) {
	
	
	// Retrieve the farthest s coordinate in the reference points
	double farthest_s = getFrenet(ref_points[0].back(), ref_points[1].back(), car_yaw, map_x, map_y)[0];
	
	// Add evenly spaced points straight down the lane line by using frenet coordinates
	int spacing = 30;
	for (int n = 1; n < 5; n++) {
		vector<double> next_points = getXY(farthest_s+spacing*n, 2+4*ref_lane, map_s, map_x, map_y);
		ref_points[0].push_back(next_points[0]);
		ref_points[1].push_back(next_points[1]);
	}
	
	// Convert to vehicle space to make sure that nothing weird happens with the math
	for (int n = 0; n < ref_points[0].size(); n++) {
		double x_ref = ref_points[0][n] - car_x;
		double y_ref = ref_points[1][n] - car_y;
		ref_points[0][n] = x_ref*cos(-deg2rad(car_yaw))-y_ref*sin(-deg2rad(car_yaw));
		ref_points[1][n] = x_ref*sin(-deg2rad(car_yaw))+y_ref*cos(-deg2rad(car_yaw));
	}
	
	// Create a spline on the reference points
	tk::spline s;
	s.set_points(ref_points[0], ref_points[1]);
	
	return s;
}

// Returns an xy trajectory given a spline and a reference velocity
vector<vector<double>> get_trajectory(tk::spline s, double ref_v, const double& car_x, const double& car_y, const double& car_yaw, const double& car_speed, const std::vector<double>& prev_path_x, const std::vector<double>& prev_path_y) {

	// Starting with the previous path for output
	vector<double> next_x_values = prev_path_x;
	vector<double> next_y_values = prev_path_y;
	
	// if we have a previous path, use it to reconstruct the history variables
	double last_x;
	double last_y;
	double planned_speed;
	if (prev_path_x.size() > 1) {
		last_x = prev_path_x.back();
		last_y = prev_path_y.back();
		planned_speed = round((sqrt((last_x - prev_path_x[prev_path_x.size()-2])*(last_x - prev_path_x[prev_path_x.size()-2]) + (last_y - prev_path_y[prev_path_y.size()-2])*(last_y - prev_path_y[prev_path_y.size()-2]))/0.02*2.23694)/0.224)*0.224;
	} else {
		last_x = car_x;
		last_y = car_y;
		planned_speed = 0.0;
	}
	
	
	
	// Fill in however many points are left blank since last update
	for (int i = 1; i <= 50-prev_path_x.size(); i++) {
		
		// Check to see if we're at the desired velocity yet
		if (planned_speed/2.23694 < ref_v) {
			planned_speed = min(planned_speed+0.224, ref_v);
		} else if (planned_speed/2.23694 > ref_v) {
			planned_speed = max(planned_speed-0.224, ref_v);
		}
		
		// Convert to meters x for this increment
		double x_increment = planned_speed*0.02/2.23694;
		
		// Copy the last xy coordinates to vehicle space
		double vehicle_last_x = last_x - car_x;
		double vehicle_last_y = last_y - car_y;
		vehicle_last_x = vehicle_last_x*cos(-deg2rad(car_yaw))-vehicle_last_y*sin(-deg2rad(car_yaw));
		
		// Increment the point
		double x_point = vehicle_last_x + x_increment;
		double y_point = s(x_point);
		
		// Convert to map space
		double ref_x = x_point*cos(deg2rad(car_yaw)) - y_point*sin(deg2rad(car_yaw)) + car_x;
		y_point = x_point*sin(deg2rad(car_yaw)) + y_point*cos(deg2rad(car_yaw)) + car_y;
		x_point = ref_x;
		
		// Append to final list
		next_x_values.push_back(x_point);
		next_y_values.push_back(y_point);
		
		// Update increment variables for next loop
		last_x = x_point;
		last_y = y_point;
		
	}
	
	return {next_x_values, next_y_values};
	
}


/* *********************************************************************
 * 
 *                          Planner Class 
 * 
 * ********************************************************************/

class Planner {
	
	public:
	
		// History variables
		int prior_plan;
		float prior_plan_discount;
	
		// Constructor
		Planner() {
			prior_plan = 5; // Start at index 5 as it's to accelerate in the center lane
			prior_plan_discount = 0.0; // Discount the effect of prior plan costs over time so as to not get stuck in a local minimum
			};
		
		// Destructor	
		~Planner() {};
		
		int choose_path(const vector<vector<vector<double>>>& possible_trajectories, const int& prev_path_size, const vector<vector<double>>& sensor_fusion, const vector<double>& map_x, const vector<double>& map_y, const vector<double>& map_s, const double& car_speed, const double& car_yaw, int& prior_plan, float& prior_plan_discount) {
			
			// Create placeholders for costs
			vector<double> costs;
			for (int i = 0; i < possible_trajectories.size(); i++) {
				costs.push_back(0.0);
			}

			/* ***********************************************************
			 *                      Cost Parameters
			 * **********************************************************/

			// Speed cost parameters
			double target_speed = 48.0;
			double speed_limit = 50.0;
			double stop_cost = 0.8;
			double speed_coef = 1.0;
			
			// Acceleration cost coef
			double accel_coef = 500.0;
			
			// Jerk cost coef
			double jerk_coef = 500.0;

			// Switching plan cost
			double prior_plan_coef = 20.0;
			double double_lane_coef = 250.0;
			
			// buffer / collusion cost
			double collusion_coef = 250.0;
			double follow_ratio = 0.225;
			
			/* ***********************************************************
			 *                     END Cost Parameters
			 * **********************************************************/
			
			// Only calculate costs for the newly created waypoints... the old waypoints will be the same for every possible plan
			double allowable_distance = 0.35;
			for (int i=0; i < possible_trajectories[0][0].size(); i++) {
				
				// TODO: Make better predictions
				// Make predictions
				double delta_t = (i+1)*0.02;
				vector<vector<double>> vehicle_predictions = sensor_fusion;
				for (int x = 0; x < vehicle_predictions.size(); x++) {
					double dist = sqrt((vehicle_predictions[x][3]*delta_t)*(vehicle_predictions[x][3]*delta_t)+(vehicle_predictions[x][4]*delta_t)*(vehicle_predictions[x][4]*delta_t));
					vehicle_predictions[x][5] = vehicle_predictions[x][5]+dist;
					}
					
					/* ***********************************************************
					 *                           COSTS
					 * **********************************************************/
					for (int x = 0; x < possible_trajectories.size(); x++) {
						
						// Determine current speed at this timestep
						double current_speed;
						if (i > 0) {
							double delta_x = possible_trajectories[x][0][i]-possible_trajectories[x][0][i-1];
							double delta_y = possible_trajectories[x][1][i]-possible_trajectories[x][1][i-1];
							current_speed = sqrt(delta_x*delta_x+delta_y*delta_y)/0.02*2.23694;
						} else {
							current_speed = car_speed;
						}
						
						// Apply speed cost
						if (current_speed < target_speed) {
							costs[x] += (stop_cost*(round(target_speed-current_speed)/target_speed))*speed_coef;
						}
						else if (current_speed > speed_limit) {
							costs[x] += 1.0*speed_coef;
						} else {
							costs[x] += (round(current_speed-target_speed)/(speed_limit-target_speed))*speed_coef;
						}
						
						// If able to calculate an acceleration, apply acceleration cost
						double prior_delta_x;
						double prior_delta_y;
						double prior_speed;
						double current_acceleration;
						if (i > 1) {
							prior_delta_x = possible_trajectories[x][0][i-1]-possible_trajectories[x][0][i-2];
							prior_delta_y = possible_trajectories[x][1][i-1]-possible_trajectories[x][1][i-2];
							prior_speed = sqrt(prior_delta_x*prior_delta_x+prior_delta_y*prior_delta_y)/0.02*2.23694;
							current_acceleration = fabs(current_speed - prior_speed);
							if (current_acceleration > 2.0) {
								costs[x] += accel_coef;
							}
						}
						
						// If able to calculate jerk, apply jerk cost
						if (i > 2) {
							double old_delta_x = possible_trajectories[x][0][i-1]-possible_trajectories[x][0][i-2];
							double old_delta_y = possible_trajectories[x][1][i-1]-possible_trajectories[x][1][i-2];
							double old_speed = sqrt(old_delta_x*old_delta_x+old_delta_y*old_delta_y)/0.02*2.23694;
							double prior_acceleration = fabs(prior_speed - old_speed);
							double current_jerk = fabs(current_acceleration - prior_acceleration);
							if (current_jerk > 2.0) {
								costs[x] += jerk_coef;
							}
						}
												
						// Apply cost for switching plans
						double discount_value = (-1/(1+exp(-0.01*prior_plan_discount))+1)*2; // the longer the car stays in the lane, the less cost to leave it
						if (!((((x == 0) | (x == 1) | (x == 2)) & ((prior_plan == 0) | (prior_plan == 1) | (prior_plan == 2))) |
						      (((x == 3) | (x == 4) | (x == 5)) & ((prior_plan == 3) | (prior_plan == 4) | (prior_plan == 5))) |
						      (((x == 6) | (x == 7) | (x == 8)) & ((prior_plan == 6) | (prior_plan == 7) | (prior_plan == 8))))) {
							costs[x] += discount_value*prior_plan_coef; // no cost for switching speeds, only for switching to a new lane
						}
						
						// Apply cost for double lane changes... pretty much just an emergency thing
						if ((((x == 0) | (x == 1) | (x == 2)) & ((prior_plan == 6) | (prior_plan == 7) | (prior_plan == 8))) |
						      (((x == 6) | (x == 7) | (x == 8)) & ((prior_plan == 0) | (prior_plan == 1) | (prior_plan == 2)))) {
							costs[x] += double_lane_coef;
						}
						
						
						// Apply cost for collusions / buffer
						bool changing_lanes = !((((x == 0) | (x == 1) | (x == 2)) & ((prior_plan == 0) | (prior_plan == 1) | (prior_plan == 2))) |
		                                (((x == 3) | (x == 4) | (x == 5)) & ((prior_plan == 3) | (prior_plan == 4) | (prior_plan == 5))) |
					                          (((x == 6) | (x == 7) | (x == 8)) & ((prior_plan == 6) | (prior_plan == 7) | (prior_plan == 8))));
						vector<double> coords = getFrenet(possible_trajectories[x][0][i], possible_trajectories[x][1][i], car_yaw, map_x, map_y);
						double current_lane = (coords[1]-2)/4;
						double backwards_final_time_to_collusion = 9999999999999.0; // Just using a big number beyond sensor range
						double forwards_final_time_to_collusion = 9999999999999.0; // Just using a big number beyond sensor range
						for (int z = 0; z < vehicle_predictions.size(); z++) {
							double backwards_time_to_collusion = 9999999999999.0; // actually this is distance now... TODO: clean up variable names
							double forwards_time_to_collusion = 9999999999999.0; // actually this is distance now... TODO: clean up variable names
							double other_vehicle_lane = (vehicle_predictions[z][6]-2)/4;
							if (abs(current_lane-other_vehicle_lane) <= 0.5) { // if overlapping lanes...
								if ((vehicle_predictions[z][5]-coords[0] <= 4.8) & (changing_lanes)) {
									backwards_time_to_collusion = sqrt((vehicle_predictions[z][5]-coords[0])*(vehicle_predictions[z][5]-coords[0]) +
																									(vehicle_predictions[z][6]-coords[1])*(vehicle_predictions[z][6]-coords[1])); // only look behind if we're changing lanes
								} 
								if (vehicle_predictions[z][5]-coords[0] >= -4.8) {
								forwards_time_to_collusion = sqrt((vehicle_predictions[z][5]-coords[0])*(vehicle_predictions[z][5]-coords[0]) +
																								(vehicle_predictions[z][6]-coords[1])*(vehicle_predictions[z][6]-coords[1]));
								} 
								if (backwards_time_to_collusion < backwards_final_time_to_collusion) {
									backwards_final_time_to_collusion = backwards_time_to_collusion;
								}
								if (forwards_time_to_collusion < forwards_final_time_to_collusion) {
									forwards_final_time_to_collusion = forwards_time_to_collusion;
								}
							}	
						}
						double effective_collusion_coef = collusion_coef;
						if (changing_lanes) {effective_collusion_coef *= 20;} // extra costs if we're changing lanes... make sure it's clear
						if (backwards_final_time_to_collusion < 9999999999999.0) {
							costs[x] += (2*(1/(1+exp(follow_ratio*backwards_final_time_to_collusion))))*effective_collusion_coef; // apply buffer cost on the closest vehicle behind with predictions
						}
						if (forwards_final_time_to_collusion < 9999999999999.0) {
							costs[x] += (2*(1/(1+exp(follow_ratio*forwards_final_time_to_collusion))))*effective_collusion_coef; // apply buffer cost on the closest vehicle ahead with predictions
						}
					}
				}
			
			// ArgMin and update the prior plan
			int chose_plan = min_element(costs.begin(), costs.end()) - costs.begin();
			if (((((chose_plan == 0) | (chose_plan == 1) | (chose_plan == 2)) & ((prior_plan == 0) | (prior_plan == 1) | (prior_plan == 2))) |
		       (((chose_plan == 3) | (chose_plan == 4) | (chose_plan == 5)) & ((prior_plan == 3) | (prior_plan == 4) | (prior_plan == 5))) |
					 (((chose_plan == 6) | (chose_plan == 7) | (chose_plan == 8)) & ((prior_plan == 6) | (prior_plan == 7) | (prior_plan == 8))))) {
				prior_plan_discount += 1;
			} else {
				prior_plan_discount = 0;
			}
			prior_plan = chose_plan;
			return prior_plan;
		}
		
		// Plan a path using features about the car, road, and other vehicles
		vector<vector<double>> plan(const double& car_x, const double& car_y, const double& car_s, const double& car_d, const double& car_yaw, const double& car_speed,
		                            const vector<double>& previous_path_x, const vector<double>& previous_path_y, const double& end_path_s, const double& end_path_d,
		                            const vector<double>& map_waypoints_x, const vector<double>& map_waypoints_y, const vector<double>& map_waypoints_s, const vector<double>& map_waypoints_dx, const vector<double>& map_waypoints_dy,
		                            const vector<vector<double>>& sensor_fusion) {
			
			// Get a pair of starter reference points near the origin of the car that every trajectory will use in the spline
			vector<vector<double>> starter_points =  get_starter_reference_points(previous_path_x, previous_path_y, car_x, car_y, car_yaw);
			
			// Generate a spline for each possible lane (note these are in car coordinates!)
			vector<tk::spline> splines;
			for (int i = 0; i < 3; i++) {
				splines.push_back(get_spline(starter_points, i, car_x, car_y, car_s, car_yaw, map_waypoints_x, map_waypoints_y, map_waypoints_s));
			}
			
			// Generate a trajectory for each direction (lane 1, 2, 3) at each speed (minimum, constant, maximum)
			vector<vector<vector<double>>> possible_trajectories;
			vector<double> speeds = {2.0, car_speed, 48.0};
			for (int t = 0; t < splines.size(); t++) {
				
				// Get the spline
				tk::spline traj = splines[t];
				
				for (int s = 0; s < speeds.size(); s++) {
					
					// Get the speed
					double spee = speeds[s];
					
					// Get the trajectory and append to the list of possibilities
					possible_trajectories.push_back(get_trajectory(traj, spee, car_x, car_y, car_yaw, car_speed, previous_path_x, previous_path_y));
					
				}
			}
			
			int chosen_plan = choose_path(possible_trajectories, previous_path_x.size(), sensor_fusion, map_waypoints_x, map_waypoints_y, map_waypoints_s, car_speed, car_yaw, prior_plan, prior_plan_discount);
			return possible_trajectories[chosen_plan];
		};
};

/* *********************************************************************
 * 
 *                             Main loop 
 * 
 * ********************************************************************/

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);
  
  // Create a planner
  Planner planner;

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&planner](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	vector<double> previous_path_x = j[1]["previous_path_x"];
          	vector<double> previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;
          	
          	// if we have more than 2 of the previous plan leftover, just drop so we can replan
          	if (previous_path_x.size() > 2) {
							previous_path_x.resize(2);
							previous_path_y.resize(2);
						}
          	
          	// define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	vector<vector<double>> plan = planner.plan(car_x, car_y, car_s, car_d, car_yaw, car_speed,
		                            previous_path_x, previous_path_y, end_path_s, end_path_d,
		                            map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy,
		                            sensor_fusion);
          	
						// Submit the trajectory
          	msgJson["next_x"] = plan[0];
          	msgJson["next_y"] = plan[1];

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































