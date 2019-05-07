#include <string.h>
#include <z3++.h>
#include <vector>
#include <map>
#include <stack>
#include <unordered_set>
#include <set>
#include <queue>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <assert.h>

using namespace std;


class EESampler {
	std::string input_file;

	struct timespec start_time;
	double solver_time = 0.0;
	int max_samples;
	double max_time;
	bool debug;
	bool random_;
	bool flip_;

	z3::context c;
	z3::optimize opt;
	std::vector<int> ind;
	std::vector<int> var;
	int samples = 0;
	int solver_calls = 0;

	std::ofstream results_file;


	///////////////////////////////////////
	// add variables
	///////////////////////////////////////
	unordered_set<int> backbone; 
	vector<unordered_set<int>> clauses;  // clauses with the form [(1,2,-3), (-1,2,3)]
	vector<unordered_set<int>> clauses_repeat;//same as clauses, but uesed for find backbone
	unordered_map<int, unordered_set<int>> clause_of_literal;   //  clauses set of a specific litearl
	unordered_map<int, unordered_set<int*>> clause_of_literal_repeat;   // pointer version used for propagation
	unordered_map<int, int> ind_to_var;
	unordered_set<string> models;  // obtained models
	vector<int*> clause_pointer_collector;
	std::unordered_set<int> varset;
	int max_var = 0;
	int run_main_count = 0;
	bool first_run_flip = true;
	int fliped_solution_num = 0, derived_solution_num = 0;
	int min_ind_count ;
	float repeat_rate = 0.0;

public:
	EESampler(std::string input, int max_samples, double max_time, bool debug, bool random_, bool flip_) : opt(c), input_file(input), max_samples(max_samples), max_time(max_time), debug(debug), random_(random_), flip_(flip_){

	}

	///////////////////////////////////////////////////////
	// reuse run() for Algorithm 3
	//////////////////////////////////////////////////////
	void run() {
		backbone.clear();
		clauses.clear();
		clauses_repeat.clear();
		clause_of_literal.clear();
		clause_of_literal_repeat.clear();
		models.clear();
		//ind_to_var.clear();
		clause_pointer_collector.clear();
		clock_gettime(CLOCK_REALTIME, &start_time);
		srand(start_time.tv_sec);
		parse_cnf();   // Line 2-5
		min_ind_count = ind.size();
		if (random_ && flip_)
			results_file.open(input_file + ".samples_random_flip");
		else if (random_)
			results_file.open(input_file + ".samples_random");
		else if (flip_)
			results_file.open(input_file + ".samples_flip");
		else
			results_file.open(input_file + ".samples_normal");

		//z3::expr_vector empty(c);
		//string tmp_str = "001000001010110100000010000101000101101001000001000101110001110010010000";
		
		int unique_count = 1;
		int repeat_count = 1;
		int control_num = -1;
		queue<int> solution_flow; 
		//Xu bug: need to excluse new back bone
		while (true) {   // Line 7
			if (debug)
				cout << "main loop " << run_main_count << endl;
			opt.push();
	/*		for (int i; i < ind.size(); i++) {
				int tmp = 2 * (tmp_str[i] - '0') - 1;
				opt.add(literal(tmp*ind[i]), 1);
			}*/

			if(solution_flow.size() > 200){
				if (solution_flow.front())
					unique_count -=1;
				else
					repeat_count -= 1;
				solution_flow.pop();
			}

			repeat_rate = float(repeat_count)/float(unique_count);
			//if (run_main_count%10==0)
			//	cout << "repeat_rate:" << repeat_rate<<endl;
			if(repeat_rate >= 0.1)
				control_num += 1;
			else if(repeat_rate <= 0.02)
				control_num -= 1;

			if (run_main_count > 6 ){
				if(control_num < 3)
					control_num = 3;
				else if(control_num > min_ind_count)
					control_num = min_ind_count;
			}


			for (int l : ind) {
				if (backbone.find(l) != backbone.end() || backbone.find(-l) != backbone.end()) continue; // Line 7 Xu continue？
				if (random_ && ((rand() % min_ind_count) >= control_num)) continue;
				if (rand() % 2) {
					opt.add(literal(l), 1);
					//tmp_str += "1";
				}
				else {
					opt.add(!literal(l), 1);
					//tmp_str += "0";
				}
				
			}
			
			if (!solve()) {
				std::cout << "Could not find a solution!\n";
				results_file << "Could not find a solution!" << endl;
				finish();
			}
			z3::model m = opt.get_model();
			string v = model_string(m, ind);
			string v_all = model_string(m, var);
			//bool equal = (tmp_str == v);
			//if (debug && equal)
			//	cout << "bingo random guess~" << endl;
			opt.pop();
			bool unique = write_solution(v);
			if (unique) {
				unique_count += 1;
				solution_flow.push(1);
				if (debug)
					cout << "find solution " << v << endl;
				if (flip_ || run_main_count < 1)
					flip(v, v_all);
				else 
					derive(v_all, v);
				print_stats(false);
			}
			else {
				repeat_count += 1;
				solution_flow.push(0);
				if (debug)
					cout << "find repeat solution in main" << v << endl;
			}
			run_main_count++;
		}

		print_stats(false);
		finish();
	}

	bool write_solution(string& v) {
		if (models.find(v) == models.end()) {
			models.insert(v);
			results_file << models.size() << " : " << v << endl;
			if (models.size() > max_samples)
				finish();
			return true;
		}
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	// Algorithm 1 with cnt
	/////////////////////////////////////////////////////////////////////////
	// int can be both positive and negative, 12 is equivalent to <12, 1>, -12 is equivalent to <12,0>
	void propagate(string v, int x) {
		if (debug)
			cout << "enter propagation" << endl;
		stack<int> B;
		
		B.push(x);
		while (!B.empty()) {
			int l = B.top();   // Line 3
			B.pop();
			for (int* p : clause_of_literal_repeat[l]) {//Line 5 Xu modified to change the value of all pointers to -1
				*p = -1;
			}
			for (int* p : clause_of_literal_repeat[-l]) {   // Line 6
				int clause_id = *p;
				if (clause_id < 0)
					continue;
				clauses_repeat[clause_id].erase(-l);
				if (clauses_repeat[clause_id].size() == 1) {  // Line 8
					*p = -1;  // Line 9
					int bb;
					for(int tmp: clauses_repeat[clause_id])
						bb = tmp;
					if (backbone.find(bb) == backbone.end()) {
						B.push(bb);  // Line 12, Fu:  we always push the variable in to B
						backbone.insert(bb);  // Line 13 to Line 16
						opt.add(literal(bb));
						if (debug)
							cout << "propagate backbone " << bb << endl;
					}
				}
			}
		}
		if (debug)
			cout << "end propagation " << endl;
	}

	///////////////////////////////////////////////////////////
	// Algorithm 2 Derivation
	//////////////////////////////////////////////////////////
	// Solutions is the parameter 'models'
	void derive(string&  v_all, string v_ind) {
		if (debug)
			cout << "enter derivation" << endl;
		queue<string> Q;
		Q.push(v_ind);
		queue<unordered_set<int>> changed_bits_queue;
		unordered_set<int> changed_bits;
		changed_bits_queue.push(changed_bits);
		unordered_map<int, int> ones_counter;
		string new_v;
		int new_derivation_number = 0;

		for (int i = 0; i < var.size(); i++) {
			int tmp = 2 * (v_all[i] - '0') - 1;
			int l = tmp*var[i];
			for (int clause : clause_of_literal[l]) {
				if (ones_counter.find(clause) == ones_counter.end())
					ones_counter.insert(make_pair(clause, 1));
				else
					ones_counter[clause] += 1;
			}
		}
		unordered_map<int, int> tmp_ones_counter;

		bool stop = false;
		while (!Q.empty() && !stop) {
			if (debug)
				cout << Q.size() << " need to be derived" << endl;
			changed_bits = changed_bits_queue.front();
			changed_bits_queue.pop();
			v_ind = Q.front();
			Q.pop();
			tmp_ones_counter = ones_counter;
			for (int bit : changed_bits) {
				int tmp = 2 * (v_ind[bit] - '0') - 1;
				int l = tmp*ind[bit];
				for (int clause : clause_of_literal[l])
					tmp_ones_counter[clause] += 1;
				for (int clause : clause_of_literal[-l])
					tmp_ones_counter[clause] -= 1;
			}

			for (int i = 0; i < ind.size(); i++) {
				bool can_derive = true;
				int tmp = 2 * (v_ind[i] - '0') - 1;
				int l = tmp*ind[i];
				if (backbone.find(l) != backbone.end())
					continue;
				for (int clause : clause_of_literal[l]) {
					if (tmp_ones_counter[clause] <= 1) {
						can_derive = false;
						break;
					}
				}
				if (can_derive) {
					new_v = v_ind;
					new_v.replace(i, 1, 1, char('0' + char(v_ind[i] == '0')));
					if (models.find(new_v) == models.end()) {
						write_solution(new_v);
						new_derivation_number += 1;
						if (debug) {
							cout << "derived " << new_v << endl;
							string new_v_all = v_all;
							if (varset.find(ind[i]) != varset.end())
								new_v_all.replace(ind_to_var[i], 1, 1, new_v[i]);
							assert(check_answer(new_v_all));
						}
						if (new_derivation_number >= 100000){
							stop = true;
							break;
						}
						Q.push(new_v);
						unordered_set<int> new_changed_bits = changed_bits;
						if (new_changed_bits.find(i) == new_changed_bits.end())
							new_changed_bits.insert(i);
						else
							new_changed_bits.erase(i);
						changed_bits_queue.push(new_changed_bits);
					}
					else if (debug)
						cout << "derived repeated solution" << endl;
				}
			}
		}

		derived_solution_num += new_derivation_number;
		if (debug)
			cout << "end derivation" << endl;
	}

	////////////////////////////////////////////
	// Algorithm 4 Flip
	///////////////////////////////////////////
	void flip(string v, string& v_all) {
		derive(v_all, v);
		if (debug)
			cout << "enter flip" << endl;	

		opt.push();
		if (flip_) {
			for (int i; i < ind.size(); i++) {
				int tmp = 2 * (v[i] - '0') - 1;
				int x = tmp * ind[i];
				opt.add(literal(x), 1);
			}
		}


		for (int i; i < ind.size(); i++) {
			int tmp = 2 * (v[i] - '0') - 1;
			int x = tmp * ind[i];
			if (backbone.find(x) != backbone.end()){
				if(first_run_flip)
					min_ind_count -= 1;
				continue;
			}
			bool solved = false;
			string new_v = v, new_v_all;
			new_v.replace(i, 1, 1, char('0' + char(v[i] == '0')));
			if (models.find(new_v) != models.end())
				continue;
			opt.push();
			opt.add(literal(-x));
			//z3::expr_vector tmp_vec(c);
			//tmp_vec.push_back(literal(-x));
			if (solved = solve()) {//give solved a value
				new_v = model_string(opt.get_model(), ind);
				new_v_all = model_string(opt.get_model(), var);
			}
			opt.pop();

			if (solved) {
				bool unique = write_solution(new_v);

				if (unique) {
					fliped_solution_num += 1;
					if (debug) {
						assert(check_answer(new_v_all));
						cout << "find flip solution at " << i << ' ' << new_v << endl;
					}
				}
				else if(debug)
					cout << "repeat flip solution at " << i << ' ' << new_v << endl;
				//Xu delete if arg_ind
				derive(new_v_all, new_v);
			}
			else if (first_run_flip) {
				if (debug)
					cout << "find backbone " << x << endl;
				min_ind_count -= 1;
				backbone.insert(x);
				opt.add(literal(x));
				propagate(v, x);
			}
		}
		opt.pop();

		if (first_run_flip) {
			for (int bb : backbone)
				opt.add(literal(bb));
			first_run_flip = false;
		}
		if (debug)
			cout << "end flip" << endl;
	}


	bool check_answer(string v_all) {
		if (debug)
			cout << "enter check_answer" << endl;

		unordered_map<int, int> ones_counter;
		for (int i = 0; i < clauses.size(); i++) {
			ones_counter[i] = 0;
		}

		for (int i = 0; i < var.size(); i++) {
			int tmp = 2 * (v_all[i] - '0') - 1;
			int l = tmp*var[i];
			for (int clause : clause_of_literal[l]) {
					ones_counter[clause] += 1;
			}
		}
		for (auto& tmp : ones_counter) {
			if (tmp.second == 0)
				return false;
		}
		return true;
	}

	void print_stats(bool simple) {
		struct timespec end;
		clock_gettime(CLOCK_REALTIME, &end);
		double elapsed = duration(&start_time, &end);
		results_file << "Samples " << models.size() << '\n';
		results_file << "Execution time " << elapsed << '\n';
		if (simple)
			return;
		results_file << "Solver time: " << solver_time << '\n';
		results_file << "Solver calls: " << solver_calls << '\n';
		results_file << "Fliped number: " << fliped_solution_num << '\n';
		results_file << "Derived number " << derived_solution_num << '\n';
		results_file << "Backbones number " << backbone.size() << '\n';
		results_file << "repeat rate " <<  repeat_rate << "\n";
		results_file << "\n";
	}

	void parse_cnf() {
		if (debug)
			cout << "enter parse" << endl;
		z3::expr_vector exp(c);
		std::ifstream f(input_file);
		if (!f.is_open()) {
			std::cout << "Error opening input file\n";
			abort();
		}
		std::unordered_set<int> indset;
		bool has_ind = false;
		std::string line;
		int counter = 0;
		int* clause_pointer = new int{ counter };
		clause_pointer_collector.push_back(clause_pointer);
		clause_of_literal.clear();
		clause_of_literal_repeat.clear();
		unordered_set<int> cls;
		int tmp;
		std::string tmps;
		while (getline(f, line)) {
			std::istringstream iss(line);
			if (line.find("c ind ") == 0) {
				iss >> tmps;
				iss >> tmps;
				while (!iss.eof()) {       //// using -i to enable the independent support
					iss >> tmp;
					if (tmp && indset.find(tmp) == indset.end()) {
						indset.insert(tmp);
						ind.push_back(tmp);
						has_ind = true;
					}
				}
			}

			else if (line[0] != 'c' && line[0] != 'p') {
				// add clauses from CNF to vector<unordered_set<int>> clauses
				if (debug && counter % 2000 == 0)
					cout << counter << " lines of cnf parsed" << endl;
				z3::expr_vector clause(c);
				while (!iss.eof()) {
					iss >> tmp;
					if (tmp != 0) {
						clause.push_back(literal(tmp));
						cls.insert(tmp);
						varset.insert(abs(tmp));//add by Xu
						if (!has_ind)
							indset.insert(abs(tmp));
						if (clause_of_literal.find(tmp) == clause_of_literal.end()) {
							unordered_set<int*> clause_ids_repeat{ clause_pointer };
							unordered_set<int> clause_ids{ counter };
							clause_of_literal_repeat.insert(std::make_pair(tmp, clause_ids_repeat));
							clause_of_literal.insert(std::make_pair(tmp, clause_ids));
						}
						else {
							clause_of_literal_repeat[tmp].insert(clause_pointer);
							clause_of_literal[tmp].insert(counter);
						}
					}
					else if (clause.size() > 0) {
						counter++;//Xu modified there may be more or less than one clause in one line
						clause_pointer = new int{ counter };
						clause_pointer_collector.push_back(clause_pointer);
						exp.push_back(mk_or(clause));
						clauses.push_back(cls);
						clauses_repeat.push_back(cls);
						clause.push_back(literal(0));
						cls.clear();
					}					
				}
			}
		}

		if (debug)
		{
			cout << counter << " lines of cnf parsed" << endl;
			cout << "file read end" << endl;
			cout << varset.size() << " vars and " << indset.size() << " inds and " << clauses.size() << " clauses found" <<endl;
		}
		f.close();

		if (!has_ind) {
			for (int lit : indset)
				ind.push_back(lit);
		}

		for (int lit : varset) {
			var.push_back(lit);
		}

		if (debug)
			cout << "begin make formula" << endl;
		z3::expr formula = mk_and(exp);
		opt.add(formula);
		//int ind_size = ind.size();
		//map<int, double>weight;
		//for (int i; i < ind_size; i++) {
		//	int _ = ind[i];
		//	int ones = clause_of_literal[_].size();
		//	int zeros = clause_of_literal[-_].size();
		//	double value = (ones + 1)*(zeros + 1);
		//	weight.insert(make_pair(_, value));
		//}
		//sort(ind.begin(), ind.end(), [&](int a, int b) { return weight[a]>=weight[b]; });

		if (debug) {
			cout << "enter ind to var process" << endl;
			for (int i = 0; i < ind.size(); i++) {
				int tmp = ind[i];
				for (int j = 0; j < var.size(); j++) {
					if (var[j] == tmp)
						ind_to_var.insert(make_pair(i, j));
				}
			}
			cout << "out ind to var process and end parse" << endl;
		}
	}

	void finish() {
		print_stats(false);
		results_file.close();
		for (int* p : clause_pointer_collector)
			delete p;
		exit(0);
	}

	bool solve() {
		struct timespec start;
		clock_gettime(CLOCK_REALTIME, &start);
		double elapsed = duration(&start_time, &start);
		if (elapsed > max_time) {
			std::cout << "Stopping: timeout\n";
			finish();
		}
		if (samples >= max_samples) {
			std::cout << "Stopping: samples\n";
			finish();
		}
		usleep(1);
		z3::check_result result = opt.check();
		struct timespec end;
		clock_gettime(CLOCK_REALTIME, &end);
		solver_time += duration(&start, &end);
		solver_calls += 1;

		return result == z3::sat;
	}


	std::string model_string(z3::model model, vector<int>& vec) {
		std::string s;

		for (int v : vec) {
			z3::func_decl decl(literal(v).decl());
			z3::expr b = model.get_const_interp(decl);
			if (b.bool_value() == Z3_L_TRUE) {
				s += "1";
			}
			else {
				s += "0";
			}
		}
		return s;
	}


	double duration(struct timespec * a, struct timespec * b) {
		return (b->tv_sec - a->tv_sec) + 1.0e-9 * (b->tv_nsec - a->tv_nsec);
	}

	z3::expr literal(int v) {
		if (v > 0)
			return c.constant(c.str_symbol(std::to_string(v).c_str()), c.bool_sort());
		else
			return !c.constant(c.str_symbol(std::to_string(-v).c_str()), c.bool_sort());
	}

};



int main(int argc, char * argv[]) {
	int max_samples = 10000000;   // 设置最大sample值
	double max_time = 7200.0;     // 设置最大计算时间
	if (argc < 2) {
		std::cout << "Argument required: input file\n";
		abort();
	}
	bool arg_samples = false;
	bool arg_time = false; 
	bool arg_debug = false;
	bool arg_random = false;
	bool arg_flip = false;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-n") == 0)
			arg_samples = true;
		else if (strcmp(argv[i], "-t") == 0)
			arg_time = true;
		else if (strcmp(argv[i], "-d") == 0)
			arg_debug = true;
		else if (strcmp(argv[i], "-r") == 0)
			arg_random = true;
		else if (strcmp(argv[i], "-f") == 0)
			arg_flip = true;
		else if (arg_samples) {
			arg_samples = false;
			max_samples = atoi(argv[i]);
		}
		else if (arg_time) {
			arg_time = false;
			max_time = atof(argv[i]);
		}
		//Xu delete one arg_ind
	}
	EESampler s(argv[argc - 1], max_samples, max_time, arg_debug, arg_random, arg_flip);
	s.run();
	return 0;
}









































