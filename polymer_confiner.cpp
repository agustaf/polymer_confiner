#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "./Random_Number_Generators/randgen.h"

//Define the confinement geometry, 0 is cylindrical tube
// 1 is square tube, 2 is sphere.
#define CONFINEMENT_GEOMETRY 0

using namespace std;
using namespace randgenbase;

#include "./Random_Number_Generators/randgen.cpp"

//This class has functions to roll monomer 
//joint angles with appropriate distributions.
class roll_angles {
	public:
	~roll_angles() {}
	roll_angles() {
		kappa = 0.0;
		neg_inverse_kappa = 0.0;
		inverse_kappa = 0.0;
		e_kappa = 0.0;
		e_neg_kappa = 0.0;
		e_kappa_diff = 0.0;
	}
	void set_kappa(const double kappa_in) {
		kappa = kappa_in;
		inverse_kappa = 1.0/kappa_in;
		neg_inverse_kappa = -1.0/kappa_in;
		e_kappa = exp(kappa);
		e_neg_kappa = exp(-kappa);
		e_kappa_diff = e_kappa - e_neg_kappa;
		return;
	}
	double roll_theta() {
		const double theta = acos(inverse_kappa* \
		  log(e_kappa-randgen()*e_kappa_diff));
		return(theta);
	}
	double roll_costheta() {
		const double costheta = inverse_kappa* \
		  log(e_kappa-randgen()*e_kappa_diff);
		return(costheta);
		/*
		if((costheta < -1.0) || (costheta > 1.0)){
			cout << "Error, rolled costheta out of range" << endl;
		}
		*/
	}
	double roll_phi() {
		return(randgen()*two_pi);
	}
	double roll_cosphi() {
		return(cos(randgen()*two_pi));
	}
	
	private:
	constexpr double two_pi = 6.28318530717959;
	double kappa;
	double neg_inverse_kappa;
	double inverse_kappa;
	double e_kappa;
	double e_neg_kappa;
	double e_kappa_diff;
};

int find_descriptor(ifstream& fp_in, const string descriptor) {
	std::string in_string = "empty";
	while (!fp_in.eof()) {
		fp_in >> in_string;
		if (in_string == descriptor) {
			return(0);
		}
	}
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	while (!fp_in.eof()) {
		fp_in >> in_string;
		if (in_string == descriptor) {
			return(0);
		}
	}
	cout << "find_descriptor, unable to find descriptor " \
	  << descriptor << endl;
	return(1);
}



int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "No input file specified." << endl;
		return(1);
	}
	constexpr int max_monomers = 1048576;
	constexpr int max_monomers_lessone = max_monomers-1;
	constexpr int quarter_max_monomers = int(double(max_monomers)/4.0);
	constexpr double pi = 3.14159265358979;
	
	/*
	const double kappa = 1.0;
	const double D = 3.5;
	const int trees_to_average = 2;
	*/
	double kappa(0.0);
	double D(0.0);
	int trees_to_average(0);
	string sweep_variable = "empty";
	int sweep_choice(-1);
	int sweep_count(0);
	double sweep_variable_change(0.0);
	double** sweep_results = 0;
	string inputfile = argv[1];
	ifstream fp_in;
	int input_success(0);
	
	fp_in.open(inputfile, ios::in);
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	input_success += find_descriptor(fp_in, "kappa=");
	fp_in >> kappa;
	if (kappa < 0.0) {
		cout << "Input error, negative kappa." << endl;
		input_success += 1;
	}
	input_success += find_descriptor(fp_in, "D=");
	fp_in >> D;
	if (D <= 0.0) {
		cout << "Input error, non-positive kappa." << endl;
		input_success += 1;
	}
	input_success += find_descriptor(fp_in, "trees_to_average=");
	fp_in >> trees_to_average;
	if (trees_to_average < 1) {
		cout << "Input error, non-positive trees_to_average." << endl;
		input_success += 1;
	}
	input_success += find_descriptor(fp_in, "sweep_variable=");
	fp_in >> sweep_variable;
	if (sweep_variable == "kappa") {
		sweep_choice = 0;
	}
	else if (sweep_variable == "D") {
		sweep_choice = 1;
	}
	else {
		sweep_choice = -1;
	}
	if (sweep_choice > -1) {
		input_success += find_descriptor(fp_in, "sweep_count=");
		fp_in >> sweep_count;
		if (sweep_count < 1) {
			cout << "Input error, non-positive sweep_count." << endl;
			input_success += 1;
		}
		input_success += find_descriptor(fp_in, "sweep_variable_change=");
		fp_in >> sweep_variable_change;
	}
	fp_in.close();
	if (input_success != 0) {
		cout << "Error, required input parameters  ";
		cout << "missing or out of bounds." << endl;
		return(1);
	}
	if (sweep_choice > -1) {
		sweep_results = new double*[sweep_count];
		for (int i=0; i<sweep_count; ++i) {
			sweep_results[i] = new double[2];
			sweep_results[i][0] = 0.0;
			sweep_results[i][1] = 0.0;
		}
	}
	
	const int branching_number = 4;
	const int branch_length = 55;
	const int polymer_length_in_branches = 15;
	const double polymer_length = \
	  double(polymer_length_in_branches)*double(branch_length) + 1.0;
	const double log_possible_survivors = log(double(trees_to_average)) + \
	  double(polymer_length_in_branches)*log(double(branching_number));
	const int possible_survivors = trees_to_average* \
	  pow(double(branching_number), double(polymer_length_in_branches));
	const double rounding_tolerance = 0.0005;
	
	if (sweep_choice < 0) {
		cout << "No sweep variable recognized, performing no sweeps." << endl;
		cout << "kappa = " << kappa << endl;
		cout << "D = " << D << endl;
	}
	else if (sweep_choice == 0) {
		cout << "Performing kappa sweep." << endl;
		cout << "kappa start value = " << kappa << endl;
		cout << "Number of kappa values being tried = " << sweep_count << endl;
		cout << "kappa change between simulations = " << \
		  sweep_variable_change << endl;
		cout << "D = " << D << endl;
	}
	else if (sweep_choice == 1) {
		cout << "Performing D sweep." << endl;
		cout << "D start value = " << D << endl;
		cout << "Number of D values being tried = " << sweep_count << endl;
		cout << "D change between simulations = " << \
		  sweep_variable_change << endl;
		cout << "kappa = " << kappa << endl;
	}
	cout << "branch_length = " << branch_length << endl;
	cout << "polymer_length_in_branches = " << \
	  polymer_length_in_branches << endl;
	cout << "polymer length total = " << \
	  (branch_length*polymer_length_in_branches) + 1 << endl;
	cout << "branching_number = " << branching_number << endl;
	#if CONFINEMENT_GEOMETRY == 1
	cout << "Confining to a square tube." << endl;
	#elif CONFINEMENT_GEOMETRY == 2
	cout << "Confining to a sphere." << endl;
	#else
	cout << "Confining to a cylindrical tube." << endl;
	#endif
	cout << "Beginning polymer confinement simulation." << endl;
	
	int list_count[2];
	list_count[0] = 0;
	list_count[1] = 0;
	double**** ends_list = 0;
	double** growth_pair = 0;
	
	//This allocates the memory for listA, listB, and growth_pair.
	ends_list = new double***[2];
	for (int i=0; i<2; ++i) {
		ends_list[i] = new double**[max_monomers];
		for (int j=0; j<max_monomers; ++j) {
			ends_list[i][j] = new double*[2];
			for (int k=0; k<2; ++k) {
				ends_list[i][j][k] = new double[3];
				for (int l=0; l<3; ++l) {
					ends_list[i][j][k][l] = 0.0;
				}
			}
		}
	}
	growth_pair = new double*[3];
	for (int i=0; i<3; ++i) {
		growth_pair[i] = new double[3];
		for (int j=0; j<3; ++j) {
			growth_pair[i][j] = 0.0;
		}
	}
	
	/*
	kappa = 1.0;
	roller.set_kappa(kappa);
	double numstore[10000];
	for(int i=0; i<10000; ++i){
		numstore[i] = roller.roll_theta();
	}
	ofstream fp_out;
	fp_out.open("roll_test.txt");
	for(int i=0; i<10000; ++i){
		fp_out << numstore[i] << endl;
	}
	fp_out.close();
	*/
	double costheta(0.0),sintheta(0.0);
	double cosphi(0.0),sinphi(0.0);
	int branch_fail(0);
	int tree_explosion(0);
	int list_index(0),store_index(0);
	if (sweep_count < 1) {
		sweep_count = 1;
		sweep_variable_change = 0.0;
	}
	autoinit_randgen();
	roll_angles roller;
	roller.set_kappa(kappa);
	double free_energy(0.0);
	double free_energy_per_length(0.0);
	double D_squared = D*D;
	double radius_squared = D_squared/4.0;
	double Rx = D/2.0;
	double Ry = D/2.0;
	int survivors(0);
	for (int sweep=0; sweep<sweep_count; ++sweep) {
		survivors = 0;
		for (int tree_index=0; tree_index < trees_to_average; ++tree_index) {
			branch_fail = 0;
			list_index = 0;
			store_index = 0;
			ends_list[0][0][0][0] = 0.0;
			ends_list[0][0][0][1] = 0.0;
			ends_list[0][0][0][2] = 0.0;
			ends_list[0][0][1][0] = 0.0;
			ends_list[0][0][1][1] = 0.0;
			ends_list[0][0][1][2] = 1.0;
			list_count[0] = 1;
			list_count[1] = 0;
			for (int i=0; i<polymer_length_in_branches; ++i) {
				int current_count = list_count[list_index];
				store_index = (list_index+1)%2;
				list_count[store_index] = 0;
				for (int j=0; j<current_count; ++j) {
					for (int k=0; k<branching_number; ++k) {
						growth_pair[0][0] = ends_list[list_index][j][0][0];
						growth_pair[0][1] = ends_list[list_index][j][0][1];
						growth_pair[0][2] = ends_list[list_index][j][0][2];
						growth_pair[1][0] = ends_list[list_index][j][1][0];
						growth_pair[1][1] = ends_list[list_index][j][1][1];
						growth_pair[1][2] = ends_list[list_index][j][1][2];

						int start_index = 2;
						int end_index = 0;
						int next_index = 1;
						branch_fail = 0;
						for (int l=0; l<branch_length; ++l) {
							start_index = (start_index+1)%3;
							end_index = (end_index+1)%3;
							next_index = (next_index+1)%3;
							costheta = growth_pair[end_index][2]- \
							  growth_pair[start_index][2];
							//cout << endl;
							sintheta = 1.0 - costheta*costheta;
							if (sintheta <= 0.0) {
								sintheta = 0.0;
								if (costheta > 0.0) {
									costheta = 1.0;
								}
								else {
									costheta = -1.0;
								}
								cosphi = 1.0;
								sinphi = 0.0;
							}
							else {
								sintheta = sqrt(sintheta);
								cosphi = (growth_pair[end_index][0]- \
								  growth_pair[start_index][0])/sintheta;
								sinphi = sqrt(1.0-cosphi*cosphi);
								if ((growth_pair[end_index][1]-growth_pair[start_index][1]) < 0.0) {
									sinphi = -sinphi;
								}
								//sinphi = (growth_pair[end_index][1]- \
								//  growth_pair[start_index][1])/sintheta;
								//cout << "sinphi: " <<  sinphi << " " << (growth_pair[end_index][1]-growth_pair[start_index][1])/sintheta << endl;
							}
							/*cout << "costheta = " << costheta << endl;
							cout << "sintheta = " << sintheta << endl;
							cout << "cosphi = " << cosphi << endl;
							cout << "sinphi = " << sinphi << endl;*/
							double newz = roller.roll_costheta();
							//cout << newz << endl;
							double newr = 1.0 - newz*newz;
							if (newr <= 0.0) {
								newr = 0.0;
								if (newz > 0.0) {
									newz = 1.0;
								}
								else {
									newz = -1.0;
								}
							}
							else {
								newr = sqrt(newr);
							}
							double newphi = roller.roll_phi();
							//double newphi = 0.0;
							double newx = newr*cos(newphi);
							double newy = (newr*newr)-(newx*newx);
							if (newy <= 0.0) {
								newy = 0.0;
							}
							else { 
								newy = sqrt(newy);
								if (newphi > pi) {
									newy = -newy;
								}
							}
							//double newy = newr*sin(newphi);
							//cout << "newy: " << newr*sin(newphi) << " " << newy << endl;;
							//cout << "new x,y,z = " << newx << ", " << newy << ", " << newz << endl;
							growth_pair[next_index][0] = \
							  (cosphi*costheta*newx - sinphi*newy + \
							  cosphi*sintheta*newz) + growth_pair[end_index][0];
							growth_pair[next_index][1] = \
							  (sinphi*costheta*newx + cosphi*newy + \
							  sinphi*sintheta*newz) + growth_pair[end_index][1];
							growth_pair[next_index][2] = \
							  (costheta*newz - sintheta*newx) + \
							  growth_pair[end_index][2];
							/*cout << "growth_pair at start_index: " << growth_pair[start_index][0] << " " << growth_pair[start_index][1] << " " << growth_pair[start_index][2] << endl;
							cout << "growth_pair at end_index: " << growth_pair[end_index][0] << " " << growth_pair[end_index][1] << " " << growth_pair[end_index][2] << endl;
							cout << "growth_pair at next_index: " << growth_pair[next_index][0] << " " << growth_pair[next_index][1] << " " << growth_pair[next_index][2] << endl;*/
							/*
							double lowlength_squared = 0.0;
							for(int ti=0; ti<3; ++ti){
								double lowtemp = growth_pair[end_index][ti]-growth_pair[start_index][ti];
								lowlength_squared += lowtemp*lowtemp;
							}
							double highlength_squared = 0.0;
							for(int ti=0; ti<3; ++ti){
								double hightemp = growth_pair[next_index][ti]-growth_pair[end_index][ti];
								highlength_squared += hightemp*hightemp;
							}
							*/
							//cout << "monomer lengths squared: " << lowlength_squared << " " << highlength_squared << endl;
							/*
							double monlength_diff = \
							  abs(sqrt(highlength_squared) - 1.0);
							if(monlength_diff > rounding_tolerance){
								cout << "error, monomer length outside tolerance: " << monlength_diff << endl;
							}
							double real_costheta = 0.0;
							for(int mi=0; mi<3; ++mi){
								double am = growth_pair[next_index][mi] - growth_pair[start_index][mi];
								real_costheta += am*am;
							}
							real_costheta = (real_costheta/2.0)-1.0;
							*/
							/*cout << "rolled costheta: " << newz << endl;
							cout << "applied costheta: " << real_costheta << endl;
							cout << "rolled phi: " << newphi << endl;*/
							/*
							double costheta_diff = abs(real_costheta-newz);
							if( costheta_diff > rounding_tolerance){
								cout << "error, rolled costheta - applied costheta outside tolerance: " << costheta_diff << endl;
								cout << "newz: " << newz << ", newphi: " << newphi << endl;
								cout << sinphi << " " << cosphi << " " << sinphi*sinphi+cosphi*cosphi << endl;
								cout << sintheta << " " << costheta << " " << sintheta*sintheta+costheta*costheta << endl;
								cout << newr << " " << newx << " " << newy << endl;
							}
							*/
							#if CONFINEMENT_GEOMETRY == 1
							if ((abs(growth_pair[next_index][0]) >= Rx) || (abs(growth_pair[next_index][1]) >= Ry)) {
								branch_fail = 1;
								break;
							}
							#elif CONFINEMENT_GEOMETRY == 2
							double r_squared = \
							  growth_pair[next_index][0]*growth_pair[next_index][0] + \
							  growth_pair[next_index][1]*growth_pair[next_index][1] + \
							  growth_pair[next_index][2]*growth_pair[next_index][2];
							if (r_squared >= radius_squared) {
								branch_fail = 1;
								break;
							}
							#else
							double r_squared = \
							  growth_pair[next_index][0]*growth_pair[next_index][0] + \
							  growth_pair[next_index][1]*growth_pair[next_index][1];
							if (r_squared >= radius_squared) {
								branch_fail = 1;
								break;
							}
							#endif
							/*
							int current_count = list_count[store_index];
							if(current_count > quarter_max_monomers){
								cout << "r_squared = " << r_squared << "  radius_squared = " << radius_squared << endl;
								cout << growth_pair[next_index][0] << " " << growth_pair[next_index][1] << " " << growth_pair[next_index][2] << endl;
							}
							*/
						}
						if (branch_fail == 0) {
							int newend = list_count[store_index];
							if (newend == max_monomers_lessone) {
								cout << "Abandoning tree, maximum surviving monomers reached." << endl;
								tree_explosion = 1;
								break;
							}
							ends_list[store_index][newend][0][0] = growth_pair[end_index][0];
							ends_list[store_index][newend][0][1] = growth_pair[end_index][1];
							ends_list[store_index][newend][0][2] = growth_pair[end_index][2];
							ends_list[store_index][newend][1][0] = growth_pair[next_index][0];
							ends_list[store_index][newend][1][1] = growth_pair[next_index][1];
							ends_list[store_index][newend][1][2] = growth_pair[next_index][2];
							++list_count[store_index];
						}
						else {
							branch_fail = 0;
						}
					}
					if (tree_explosion == 1) {
						break;
					}
				}
				if (tree_explosion == 1) {
					break;
				}
				else {
					list_index = (list_index+1)%2;
				}
			}
			if (tree_explosion == 1) {
				--tree_index;
				tree_explosion = 0;
				continue;
			}
			else {
				cout << "Sweep = " << sweep << ", Tree = " << tree_index << \
				  ", Tree survivors = " << list_count[store_index] << endl;
				survivors += list_count[store_index];
			}
		}
		free_energy = -log(double(survivors))+log_possible_survivors;
		free_energy_per_length = free_energy/polymer_length;
		if (sweep_choice > -1) {
			if (sweep_choice == 0) {
				sweep_results[sweep][0] = kappa;
				sweep_results[sweep][1] = free_energy_per_length;
				kappa += sweep_variable_change;
				if (kappa < 0.0) {
					cout << "Error, negative kappa in sweep, setting kappa to zero." << endl;
					kappa = 0.0;
				}
				roller.set_kappa(kappa);
			}
			else if (sweep_choice == 1) {
				sweep_results[sweep][0] = D;
				sweep_results[sweep][1] = free_energy_per_length;
				D += sweep_variable_change;
				if (D <= 0.0) {
					cout << "Error, non-positive D value in sweep, setting D to one." << endl;
					D = 1.0;
				}
				D_squared = D*D;
				radius_squared = D_squared/4.0;
				Rx = D/2.0;
				Ry = D/2.0;
			}
		}
	}
	if (sweep_choice == -1) {
		cout << endl;
		cout << "Survivors: " << survivors << endl;
		cout << "Possible Survivors: " << possible_survivors << endl;
		cout << "Free energy of confinement: " << free_energy << endl;
		cout << "Free energy per length: " << free_energy_per_length << endl;
		cout << endl;
	}
	else {
		ofstream fp_out;
		fp_out.open("outputfile.txt", ios::out);
		fp_out.clear();
		fp_out.seekp(0, ios::beg);
		#if CONFINEMENT_GEOMETRY == 1
		fp_out << "Square tube confinement simulation." << endl;
		#elif CONFINEMENT_GEOMETRY == 2
		fp_out << "Sphere confinement simulation." << endl;
		#else
		fp_out << "Cylindrical tube confinement simulation." << endl;
		#endif
		fp_out << "branch_length = " << branch_length << endl;
		fp_out << "polymer_length_in_branches = " << \
		  polymer_length_in_branches << endl;
		fp_out << "polymer length total = " << branch_length*polymer_length_in_branches+1 << endl;
		fp_out << "branching_number = " << branching_number << endl;
		if (sweep_choice == 0) {
			fp_out << "D = " << D << endl;
			fp_out << "kappa ";
		}
		else if (sweep_choice == 1) {
			fp_out << "kappa = " << kappa << endl;
			fp_out << "D ";
		}
		fp_out << "free_energy_per_length" << endl;
		for (int i=0; i<sweep_count; ++i) {
			fp_out << sweep_results[i][0] << " " << sweep_results[i][1] << endl;
		}
		fp_out.close();
	}

	//This deallocates the memory held by array structures.
	if (ends_list) {
		for (int i=0; i<2; ++i) {
			for (int j=0; j<max_monomers; ++j) {
				for (int k=0; k<2; ++k) {
					delete[] ends_list[i][j][k];
				}
				delete[] ends_list[i][j];
			}
			delete[] ends_list[i];
		}
		delete[] ends_list;
		ends_list = nullptr;
	}
	if (growth_pair) {
		for (int i=0; i<3; ++i) {
			delete[] growth_pair[i];
		}
		delete[] growth_pair;
		growth_pair = nullptr;
	}
	if (sweep_results) {
		for (int i=0; i<sweep_count; ++i) {
			delete[] sweep_results[i];
		}
		delete[] sweep_results;
		sweep_results = nullptr;
	}
	return(0);
}
