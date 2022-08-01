#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<cmath>
#include<thread>
#include<random>
#include <set>
#include<math.h>
#include<ctime>
#include <mutex>
using namespace std;
bool NextSet(vector<vector<double>>& _points);
void foo();
class experement {
public:

	experement(int id,size_t _size=1) {
		this->ID = id;
		Setnum_of_characteristics();
		copters_size = _size;
		vector<pair< int, string>>copters(copters_size);
		this->copters = copters;
		this->points_by_claster.resize(this->copters.size());
	}
	~experement() {
		copters.clear();
	}
	void Setnum_of_characteristics() {
		num_of_characteristics = 5;  //random
	}
	void randomize_points() {

	}
	void SetUpcopters(int _num) {
		int  _num_of_characteristics = 5;
		string characteristics = "";
		vector<pair< int, string>>_copters;
		vector<int>alldroncateg(_num_of_characteristics, 0);
		time_t td;
		td = time(NULL);
		srand(time_t(ctime(&td)) + this->ID);     // random
		for (int i = 0; i < _num;i++) {
			characteristics = Setcharacteristics(3,5);
			_copters.push_back(make_pair(i+1, characteristics));
		}
		for (int i = 0; i < _copters.size(); i++) {
			for (int j = 0; j < _copters[i].second.length(); j++) {  //  числа каждой категории  у дронов
				alldroncateg[_copters[i].second[j] - '0' - 1]++;
			}
		}
		for (int i = 0; i < alldroncateg.size(); i++) {
			if (alldroncateg[i]<2) {
				for (int j = 0; j < _copters.size(); j++) {  
					bool isthere = false;
					for (int l = 0; l < _copters[j].second.length();l++) {
						if ((_copters[j].second[l]-'0')==i+1) {
							isthere = true;
							break;
						}				
						else {
							
							isthere = false;
						}
					}
					
					if (!isthere) {
						_copters[j].second += to_string(i + 1);
						_copters[j].second = sortstring(_copters[j].second);
						break;
					}
				}
			}
			
		}

		Setcopters(_copters);
		
	}
	void SetUpField_points(int _num) {
		time_t td;
		td = time(NULL);
		srand(time_t(ctime(&td)) * this->ID);       // random
		vector<pair<vector<double>, string>> _field_points(_num);
		for (int i = 0; i < _num;i++) {
			for (int j = 0; j < 3;j++) {
				float r = (float)Random(0.0, 100.0);
				//r = floorf(r*10000)/10000;
				
				_field_points[i].first.push_back((double)r);
				_field_points[i].second= Setcharacteristics(1,5);
			}
		}
		Setfield_points(_field_points);
	}
	void Setcopters(vector<pair< int, string>> _copters) {
		copters = _copters;
		copters_size = _copters.size();
		this->points_by_claster.resize(this->copters.size());
		vector<double> _vec = {0.0, 0.0, 0.0 };
		
		this->copter_position.resize(this->copters.size(), make_pair(0, _vec));
	}
	void Setfield_points(vector<pair<vector<double>,string>> _field_points_xyz) {
		field_points_xyz = _field_points_xyz;
		field_points_size = _field_points_xyz.size();
	}
	void SetClasters() {
		int N = (double)this->copters_size;
		int Side = this->GetFieldxyz()[0];
		double V = Side*Side*Side;
		if (N%2==0) {
			int correction = 0.0;
			correction = Side % N;
			int sidex = Side / (N / 2);
			vector<int>_vec;
			
			for (int i = 0; i < N / 2-1;i++) {
				_vec.push_back(sidex*(i+1));
			}
			_vec.push_back(Side);
			for (int i =0; i < N / 2 - 1; i++) {
				_vec.push_back(sidex*N/2 -sidex* (i + 1));
			}
			_vec.push_back(0);
			claster_point_xyz = _vec;
		}
		else if (N%2!=0) {
			int correction = 0.0;
			correction = Side % N;
			int sidex = Side*2 /N ;
			vector<int>_vec;

			for (int i = 0; i < N / 2; i++) {
				_vec.push_back(sidex * (i + 1));
			}
			//_vec.push_back(Side);
			for (int i = 0; i < N / 2 ; i++) {
				_vec.push_back(sidex * (N / 2) - sidex * (i  ));
			}
			_vec.push_back(0);
			claster_point_xyz = _vec;
		}
	}
	string sortstring(string str) {
		vector<int>_vec;
		for (int i = 0; i < str.length(); i++) {
			_vec.push_back(str[i]-'0');
		}
		sort(_vec.begin(),_vec.end());
		string _str="";
		for (int i = 0; i < _vec.size();i++) {
			_str +=to_string( _vec[i]);
		}
		return _str;
	}
	void Sort_points_by_clasters() {
		int N = (double)this->copters_size;
		if (N%2==0) {
			
				for (int i = 0; i < field_points_xyz.size(); i++) {
					vector<vector<int>>prev = { {0, 0, 0}, { 0,0,fieldsizez }, { 0,fieldsizey,0 }};
					bool done = false;
					for (int j = 0; j < (claster_point_xyz.size())/2; j++) {
						if (field_points_xyz[i].first[0]>=0 && field_points_xyz[i].first[0] <= fieldsizex && field_points_xyz[i].first[1] >= 0 && field_points_xyz[i].first[1] <= fieldsizey && field_points_xyz[i].first[2] >= 0 && field_points_xyz[i].first[2] <= fieldsizey) {
							vector<int>nullvec = { 0,0,0 };
							vector<int>a = { claster_point_xyz[j],fieldsizey,fieldsizez };
							vector<int>b= { claster_point_xyz[j],fieldsizey,0 };
							vector<float>x= { (float)field_points_xyz[i].first[0],(float)field_points_xyz[i].first[1],(float)field_points_xyz[i].first[2] };
							if (plat_turn(nullvec,a,b,x)>0 && plat_turn(prev[0], prev[1], prev[2],x)<=0) {
							//	cout <<"claster"<<" "<<j<<" "<< field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2]<<endl;
								points_by_claster[j].push_back(field_points_xyz[i]);
								done = true;
								break;
							}
							else {
								done = false;
							}
							prev = { nullvec ,a,b };
						}
					}
					if (done) {
						continue;
					}
					/*prev = {};
					prev = { {0,0,0} ,{fieldsizex,fieldsizey,0} ,{fieldsizex,fieldsizey,fieldsizez} };*/
					for (int j = (claster_point_xyz.size()) / 2 ; j < claster_point_xyz.size(); j++) {
						if (field_points_xyz[i].first[0] >= 0 && field_points_xyz[i].first[0] <= fieldsizex && field_points_xyz[i].first[1] >= 0 && field_points_xyz[i].first[1] <= fieldsizey && field_points_xyz[i].first[2] >= 0 && field_points_xyz[i].first[2] <= fieldsizey) {
							vector<int>nullvec = { 0,0,0 };
							vector<int>a = { fieldsizex, claster_point_xyz[j],fieldsizez };
							vector<int>b = { fieldsizex, claster_point_xyz[j],0 };
							vector<float>x = { (float)field_points_xyz[i].first[0],(float)field_points_xyz[i].first[1],(float)field_points_xyz[i].first[2] };
							//cout << " " << field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2]<<" "<< plat_turn(nullvec, a, b, x) <<" "<<" "<< plat_turn(prev[0], prev[1], prev[2], x)<<" "<< claster_point_xyz[j]<<endl;
							if (plat_turn(nullvec, a, b, x) > 0 && plat_turn(prev[0], prev[1], prev[2], x) <= 0) {
							//	cout << "claster" << " " << j << " " << field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2] << endl;
								points_by_claster[j].push_back(field_points_xyz[i]);
								break;
							}
							else {
								//cout << "no";
							}
							prev = { nullvec ,a,b };
						}
					}
				}



			
		}
		else if (N%2!=0) {

			for (int i = 0; i < field_points_xyz.size(); i++) {
				vector<vector<int>>prev = { {0, 0, 0}, { 0,0,fieldsizez }, { 0,fieldsizey,0 } };
				bool done = false;
				for (int j = 0; j < (claster_point_xyz.size()) / 2; j++) {
					if (field_points_xyz[i].first[0] >= 0 && field_points_xyz[i].first[0] <= fieldsizex && field_points_xyz[i].first[1] >= 0 && field_points_xyz[i].first[1] <= fieldsizey && field_points_xyz[i].first[2] >= 0 && field_points_xyz[i].first[2] <= fieldsizey) {
						vector<int>nullvec = { 0,0,0 };
						vector<int>a = { claster_point_xyz[j],fieldsizey,fieldsizez };
						vector<int>b = { claster_point_xyz[j],fieldsizey,0 };
						vector<float>x = { (float)field_points_xyz[i].first[0],(float)field_points_xyz[i].first[1],(float)field_points_xyz[i].first[2] };
						if (plat_turn(nullvec, a, b, x) > 0 && plat_turn(prev[0], prev[1], prev[2], x) <= 0) {
							//cout << "claster" << " " << j << " " << field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2] << endl;
							points_by_claster[j].push_back(field_points_xyz[i]);
							done = true;
							break;
						}
						else {
							done = false;
						}
						prev = { nullvec ,a,b };
					}
				}
				if (done) {
					continue;
				}
				/*prev = {};
				prev = { {0,0,0} ,{fieldsizex,fieldsizey,0} ,{fieldsizex,fieldsizey,fieldsizez} };*/
				for (int j = (claster_point_xyz.size()) / 2; j < claster_point_xyz.size(); j++) {
					if (field_points_xyz[i].first[0] >= 0 && field_points_xyz[i].first[0] <= fieldsizex && field_points_xyz[i].first[1] >= 0 && field_points_xyz[i].first[1] <= fieldsizey && field_points_xyz[i].first[2] >= 0 && field_points_xyz[i].first[2] <= fieldsizey) {
						vector<int>nullvec = { 0,0,0 };
						vector<int>a = { fieldsizex, claster_point_xyz[j],fieldsizez };
						vector<int>b = { fieldsizex, claster_point_xyz[j],0 };
						vector<float>x = { (float)field_points_xyz[i].first[0],(float)field_points_xyz[i].first[1],(float)field_points_xyz[i].first[2] };
						//cout << " " << field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2]<<" "<< plat_turn(nullvec, a, b, x) <<" "<<" "<< plat_turn(prev[0], prev[1], prev[2], x)<<" "<< claster_point_xyz[j]<<endl;
						if (plat_turn(nullvec, a, b, x) > 0 && plat_turn(prev[0], prev[1], prev[2], x) <= 0) {
							//cout << "claster" << " " << j << " " << field_points_xyz[i].first[0] << " " << field_points_xyz[i].first[1] << " " << field_points_xyz[i].first[2] << endl;
							points_by_claster[j].push_back(field_points_xyz[i]);
							break;
						}
						else {
							//cout << "no";
						}
						prev = { nullvec ,a,b };
					}
				}
			}
		}
	}
	void Setcopter_position(int index,vector<double> _position) {
		copter_position[index].second = _position;
	}
	void Setcopter_pathlength(int index, double _pathlength) {
		copter_position[index].first += _pathlength;
	}
	void Dellcopter_pathlength_position(int index) {
		copter_position[index] = {};
	}
	vector<int> GetFieldxyz() {
		vector<int> arr(3);
		arr[0] = fieldsizex;
		arr[1] = fieldsizey;
		arr[2] = fieldsizez;
		return arr;
	}
	vector<pair<vector<double>, string>> Getfield_points_xyz() {
		return field_points_xyz;
	}
	double Random(double _min,double _max) {
		return (double)(rand())/RAND_MAX*(_max-_min)+_min;
	}
	float plat_turn(vector<int> a, vector<int> b, vector<int> c , vector<float> x) {
		float a21 = b[0] - a[0],
			a22 = b[1] - a[1],
			a23 = b[2] - a[2],
			a31 = c[0] - a[0],
			a32 = c[1] - a[1],
			a33 = c[2] - a[2];
		float A = { a22 * a33 - a23 * a32 },
			B = { a23 * a31 - a21 * a33 },
			C = { a21 * a32 - a22 * a31 };
		float D = -a[0] * A - a[1] * B - a[2] * C;
		return A * x[0] + B * x[1] + C * x[2] + D;
	}
	void DellPointAtSortedclasters(int clasterindex,int pointindex) {
		points_by_claster[clasterindex][pointindex] = {};
	}
	vector<vector<pair<vector<double>, string>>> GetSortedclasters() {
		
		return points_by_claster;
	}
	vector<pair< int, string>> Getcopters() {
		return copters;
	}
	double GetlengthOfway() {
		return lengthOfway;
	}
	double Getcopterpathlength(int index) {
		return copter_position[index].first;
	}
	vector<double> Getcopterposition(int index) {
		return copter_position[index].second;
	}
	vector<pair<double, vector<double>>> Getcopter_position() {
		return copter_position;
	}
private:
	int ID;
	int copters_size;
	int field_points_size;
	int num_of_characteristics;
	int fieldsizex = 100.0;
	int fieldsizey = 100.0;
	int fieldsizez = 100.0;
	double lengthOfway = 0.0;
	vector<pair<vector<double>, string>> field_points_xyz;
	vector<pair< int, string>>copters;
	vector<int>claster_point_xyz;
	vector<vector<pair<vector<double>, string>>>points_by_claster;
	vector<pair<double,vector<double>>>copter_position;
	string Setcharacteristics(int _num_of_characteristics,int _interval) {
		string characteristics = "";
		set<string>set_of_characteristics;
		
		for (int i = 0; i < 1+ rand()% (_num_of_characteristics-1 +1);i++) {
			set_of_characteristics.insert(to_string(1+rand() % (_interval -1+1) ));
		}
		set<string>::iterator it1;
		for ( it1 = set_of_characteristics.begin(); it1 != set_of_characteristics.end();it1++) {
			characteristics += *it1;
		}
		return characteristics;
	}
};

double SolveDistance(vector<double>_point1, vector<double>_point2) {
	double _length = 0.0;
	_length = sqrt( (_point1[0]- _point2[0])* (_point1[0] - _point2[0])+ (_point1[1] - _point2[1])* (_point1[1] - _point2[1])+ (_point1[2] - _point2[2])*(_point1[2] - _point2[2])) ;
	return _length;
}
double SolveDistance(vector<vector<double>> points) {
	double _length = 0.0;
	for (int i = 1; i < points.size();i++) {
		_length+= SolveDistance(points[i-1],points[i]);
	}
	return _length;
}
bool SolveDistanceByOtzig(vector<vector<double>> &points) {
	if (points.size()<=1) {
		return true;
	}
	double Tmax = 100.0;
	double Tmin = 0.0000000001;
	int T = Tmax;
	int iteration = 1;
	double length = 0.0;
	//vector<int> perestanovki;
	//for (int i = 0; i < points.size();i++) {
	//	perestanovki.push_back(i+1);
	//}
	length = SolveDistance(points);
	//cout << "pointslength old: " << length << endl;
	vector<vector<double>> fakepoints(points);
	while (NextSet(fakepoints) && iteration<=1000) {
		if (T > Tmin) {
			double savelength = length;
			double newlength = SolveDistance(fakepoints);
				if (newlength<length) {
					points = fakepoints;
					length = newlength;
				}
				else {
					double p = exp((-iteration) / 100 / (newlength) * (length));
					if (p > 1 + rand() % (99)) {

						points = fakepoints;

						length = newlength;
					}
					else{
						fakepoints = points;
					}
				}
					if (length >= newlength) {
						T = T * (length / savelength);
					}
				iteration = iteration + 1;
		}
	}
//cout << "pointslength new: " <<length <<endl;
	return true;
}
bool NextSet(vector<vector<double>> &_points) {
	int r1 = 0 +rand()%(_points.size()-1 - 0+1);
	int r2 = 0 + rand()%(_points.size()-1 - 0+1);

	vector<double>  num ( _points[r1]);
		_points[r1] = _points[r2];
		_points[r2] = num;
		return true;
}
double countscore(string str,vector<double>_vec) {
	double _score = 0.0;
	for (int i = 0; i < str.length();i++) {
		_score += _vec[str[i]-'0'-1];
	}
	return _score;
}
bool instring(string _s,string _str) {
	
	if (_s.length()<=_str.length()) {
		for (int i = 0; i < _str.length(); i++) {
			for (int j = 0; j < _s.length();j++) {
				if (_s[j]==_str[i]) {
					return true;
				}
			}
		}
	}
	else {
		for (int i = 0; i < _s.length(); i++) {
			for (int j = 0; j < _str.length(); j++) {
				if (_s[i] == _str[j]) {
					return true;
				}
			}
		}
	}
	return false;
}
int findmin(vector<pair<double, vector<double>>>_vec, vector<int> _bad) {
	int index = -1;
	double min = 1000000.0;
	for (int i = 0; i < _vec.size();i++) {
		
			if (_bad[i]==1) {
				continue;
			}
		
		if (min>_vec[i].first) {
			min = _vec[i].first;
			index = i;
			
		}
	}
	
	return index;
}
mutex thlock;
int numberthreads=0;
bool cmp(pair<int, pair<double, vector<double>>>p1, pair<int, pair<double, vector<double>>>p2) {
	return p1.second.first < p2.second.first;
}
bool compar(int p1,int p2) {
	return p1 > p2;
}
int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	unsigned int start_time = clock(); // начальное время
  

	for (int exi = 0; exi < 1000; exi++) {
		experement ex(exi);
		ex.SetUpcopters(20);
		ex.SetUpField_points(100);
		ex.SetClasters();
		ex.Sort_points_by_clasters();
		/*for (int j = 0; j < ex.GetSortedclasters().size();j++) {
			SolveDistanceByOtzig(ex.GetSortedclasters()[j]);
			cout << "done";
		}*/
		vector<int>allcateg(5, 0);
		vector<vector<double>>score_by_claster;

		/*cout << "copters" << endl;
		for (int tt = 0; tt < ex.Getcopters().size(); tt++) {
			cout << " " << ex.Getcopters()[tt].second << " ";
		}*/
		//cout << endl;
		for (int j = 0; j < ex.GetSortedclasters().size(); j++) {
			string categ = "";

			for (int tt = 0; tt < ex.GetSortedclasters()[j].size(); tt++) {  // все категории по 1 кластеру в 1 строку 
				//cout << " "<< ex.GetSortedclasters()[j][tt].second <<" ";
				categ += ex.GetSortedclasters()[j][tt].second;
			}
			//cout << endl;
			vector<int>_vec(5, 0);
			for (int _i = 0; _i < categ.size(); _i++) {
				allcateg[(categ[_i]) - '0' - 1]++; // подсчет общего числа категорий
				_vec[(categ[_i]) - '0' - 1]++;   // подсчет частных категорий
			}
			/*for (int _i = 0; _i < _vec.size(); _i++) {
				cout<<" "<<_i+1<<" "<<_vec[_i]<<" ";
			}*/
			//cout << endl;
		}
		for (int j = 0; j < ex.GetSortedclasters().size(); j++) {
			vector<double>score;
			vector<int>_vec(5, 0);
			string categ = "";

			for (int tt = 0; tt < ex.GetSortedclasters()[j].size(); tt++) {  // все категории по 1 кластеру в 1 строку 
				categ += ex.GetSortedclasters()[j][tt].second;
			}
			for (int _i = 0; _i < categ.size(); _i++) {
				_vec[(categ[_i]) - '0' - 1]++;   // подсчет частных категорий
			}
			for (int _i = 0; _i < _vec.size(); _i++) {
				if (allcateg[_i] > 0)
				{
					//cout << " " << _i + 1 << " " << (double)((double)_vec[_i] / (double)allcateg[_i]) << " ";  // подсчет отношения количества текущей категории в кластере к общему числу
					score.push_back((double)((double)_vec[_i] / (double)allcateg[_i]));
				}
				else {
					//cout << " " << _i + 1 << " " << 0 << " ";
				}
			}
			//cout << endl;
			score_by_claster.push_back(score);
		}

		vector<vector<double>>drone_score;
		for (int j = 0; j < ex.Getcopters().size(); j++) {
			double _score = 0.0;
			vector<double>_vec;
			for (int tt = 0; tt < score_by_claster.size(); tt++) {

				_score = countscore(ex.Getcopters()[j].second, score_by_claster[tt]);
				_vec.push_back(_score);
			}
			drone_score.push_back(_vec);
		}
		vector<pair<int, int>>drone_to_claster;
		double _max = 0;
		int drone_index = -1;
		int claster_index = -1;
		for (int _i = 0; _i < ex.Getcopters().size(); _i++) {
			drone_index = -1;
			claster_index = -1;
			_max = 0;
			for (int j = 0; j < drone_score.size(); j++) {

				for (int tt = 0; tt < drone_score[j].size(); tt++) {
					if ((_max < drone_score[j][tt] || (_max == drone_score[j][tt] && (drone_score[j].size() / 2 - claster_index) > (tt - drone_score[j].size() / 2))) && drone_score[j][tt] != 0.0) {
						bool yes = true;
						for (int l = 0; l < drone_to_claster.size(); l++) {
							if (drone_to_claster[l].second == tt) {

								yes = false;
								break;
							}
						}
						if (yes) {
							_max = drone_score[j][tt];
							drone_index = j;
							claster_index = tt;
						}
					}
				}
			}
			if (drone_index != -1 && claster_index != -1) {
				drone_to_claster.push_back(make_pair(drone_index, claster_index));
				drone_score[drone_index] = {};
			}
			else {
				break;
			}  // разбили дронов по первым кластерам .
		}
		/*for (int j = 0; j < drone_to_claster.size();j++) {
			cout << " "<< drone_to_claster[j].first <<" "<< drone_to_claster[j].second<< " " << endl;
		}*/
		//cout << endl;
		int _num_of_characteristics = 5;
		vector<int>alldroncateg(_num_of_characteristics, 0);
		for (int _i = 0; _i < ex.Getcopters().size(); _i++) {
			for (int j = 0; j < ex.Getcopters()[_i].second.length(); j++) {  //  числа каждой категории  у дронов
				alldroncateg[ex.Getcopters()[_i].second[j] - '0' - 1]++;
			}
		}
		vector<vector<vector<double>>> points_for_drone;
		int minus = 0;
		for (int j = 0; j < drone_to_claster.size(); j++) {
			vector<vector<double>> _vec;
			for (int tt = 0; tt < ex.GetSortedclasters()[drone_to_claster[j].second].size(); tt++) {

				if (instring(ex.GetSortedclasters()[drone_to_claster[j].second][tt].second, ex.Getcopters()[drone_to_claster[j].first].second)) {
					_vec.push_back(ex.GetSortedclasters()[drone_to_claster[j].second][tt].first);
					ex.DellPointAtSortedclasters(drone_to_claster[j].second, tt);
					minus++;
				}

			}
			points_for_drone.push_back(_vec);
		}
		for (int j = 0; j < ex.Getcopters().size() && j < points_for_drone.size(); j++) {
			SolveDistanceByOtzig(points_for_drone[j]);
		}
		for (int j = 0; j < ex.Getcopters().size() && j < points_for_drone.size(); j++) {
			ex.Setcopter_pathlength(j, SolveDistance({ 0,0,0 }, points_for_drone[j][0]));
			ex.Setcopter_position(j, points_for_drone[j][points_for_drone[j].size() - 1]);
			ex.Setcopter_pathlength(j, SolveDistance(points_for_drone[j]));
		}


		vector<vector<pair<vector<double>, string>>>_lastPoints(5);
		vector<int>countOf_lastPoints(5, 0);
		for (int j = 0; j < ex.GetSortedclasters().size(); j++) {
			for (int tt = 0; tt < ex.GetSortedclasters()[j].size(); tt++) {
				if (ex.GetSortedclasters()[j][tt].first.size() != 0) {
					countOf_lastPoints[stoi(ex.GetSortedclasters()[j][tt].second) - 1]++;
					_lastPoints[stoi(ex.GetSortedclasters()[j][tt].second) - 1].push_back(ex.GetSortedclasters()[j][tt]);
				}
			}
		}
		//vector<vector<pair<vector<double>,string>>> points_for_drone_last(ex.Getcopters().size());
		////vector<pair<double, vector<double>>>copter_position;
		//points_for_drone = {};
		//vector<int>badindex(ex.Getcopter_position().size(),0);
		//for (int j = 0;j< _lastPoints.size(); j++) {
		//	for (int tt = 0; tt < ex.Getcopter_position().size();tt++) {
		//		int _index = findmin(ex.Getcopter_position(), badindex);
		//		if (_index!=-1 && instring(ex.Getcopters()[_index].second,to_string(j+1))) {
		//			points_for_drone_last[j].push_back(_lastPoints[j][_index]);
		//			badindex[_index] = 1;
		//		}
		//	}
		//}
		vector<pair<int,pair<double, vector<double>>>>copter_position_for_sort;
		vector<vector<pair<int, pair<double, vector<double>>>>>copter_position_of_one_categ(5);
		for (int j = 0; j < ex.Getcopter_position().size(); j++) {
			copter_position_for_sort.push_back(make_pair(j,ex.Getcopter_position()[j]));
		}
		sort(copter_position_for_sort.begin(), copter_position_for_sort.end(),cmp);
		for (int tt = 0; tt < _lastPoints.size(); tt++) {
			if (_lastPoints[tt].size()>0) {
				for (int j = 0; j < copter_position_for_sort.size(); j++) {
					if (instring(ex.Getcopters()[copter_position_for_sort[j].first].second,to_string(tt+1))) {
						copter_position_of_one_categ[tt].push_back(copter_position_for_sort[j]);
					}

				}
			}
		}
		//for (int tt = 0; tt < copter_position_of_one_categ.size(); tt++) {
		//	double summ = 0.0;
		//	int T = _lastPoints[tt].size();
		//	int TT = T;
		//	for (int t = 0; t < copter_position_of_one_categ[tt].size(); t++) {
		//		summ += copter_position_of_one_categ[tt][t].second.first;
		//	}
		//	vector<int>tasksfor;
		//	int count = summ / T;
		//	int coun = T;
		//	for (int t = 0; t < copter_position_of_one_categ[tt].size(); t++) {
		//		
		//	//	int count = (1-(double)(copter_position_of_one_categ[tt][t].second.first/ summ))*T;
		//		//count %= TT;
		//		//cout << "";
		//		//TT -= count;
		//		int g = ceil(int(copter_position_of_one_categ[tt][t].second.first / count)); 
		//		coun -= g;
		//		tasksfor.push_back(g);
		//	}
		//	tasksfor[tasksfor.size() - 1] += coun;
		//	sort(tasksfor.begin(), tasksfor.end(),compar);
		//	for (int j = 0; j < copter_position_of_one_categ[tt].size();j++) {
		//		for (int jj = 0; jj < _lastPoints[tt].size();jj++) {

		//		}
		//	}
		//}
		double summ = 0.0;
		double maxx = 0;
		int minn = 1000000;
		for (int j = 0; j < ex.Getcopters().size(); j++) {
			summ += ex.Getcopterpathlength(j);
			if (maxx < ex.Getcopterpathlength(j)) {
				maxx = ex.Getcopterpathlength(j);
			}
			if (minn > ex.Getcopterpathlength(j)) {
				minn = ex.Getcopterpathlength(j);
			}
		}
		int kolvo = 0;
		for (int j = 0; j < _lastPoints.size(); j++) {
			kolvo += _lastPoints[j].size();
		}
		/*vector<vector<vector<double>>> points_for_drone_last(ex.Getcopters().size());
			pair<double,vector<double>>_copter_position;
			double _min = 100000000;
			while (_lastPoints.size()!=0) {
				int _index = -1;
				for (int j = 0; j < ex.Getcopter_position().size(); j++) {
					if (ex.Getcopter_position()[j].first < _min && ex.Getcopter_position()[j].first>=0) {
						_min = ex.Getcopter_position()[j].first;
						_index = j;
						_copter_position = ex.Getcopter_position()[j];
					}

				}
				ex.Setcopter_pathlength(_index,-1.0);
				vector<vector<double>> _vec;
				for (int t = 0; t < ex.Getcopters()[_index].second.length(); t++) {
					for (int tt = 0; tt < _lastPoints[(ex.Getcopters()[_index].second[t]) - '0' - 1].size(); tt++) {
						cout << " " << (ex.Getcopters()[_index].second[t]) - '0' - 1 << endl;
						_vec.push_back(_lastPoints[(ex.Getcopters()[_index].second[t]) - '0' - 1][tt].first);
					}
					_lastPoints[(ex.Getcopters()[_index].second[t]) - '0' - 1] = {};
				}
				points_for_drone_last[_index] = _vec;
			}*/
		
		//cout  << " " << summ + summ * kolvo / (ex.Getfield_points_xyz().size()) << endl;
		//cout <<summ/double((double)1.0-(double)kolvo / (ex.Getfield_points_xyz().size())) << endl;
		cout << maxx/ double((double)1.0 - (double)kolvo / (ex.Getfield_points_xyz().size()));
	}
	


	unsigned int end_time = clock(); // конечное время
	cout << end_time - start_time; // искомое время
	return 0;
}