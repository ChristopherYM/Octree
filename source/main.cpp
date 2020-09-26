#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

template<typename T>
struct point{
	T x, y, z;
	point(T a, T b, T c) : x(a), y(b), z(c) {}
	bool operator ==(point<T> &p){
		return (x == p.x && y == p.y && z == p.z);
	}
	bool operator !=(point<T> &p){
		return (x != p.x || y != p.y || z != p.z);
	}
};


template<typename T>
struct node{
	point<T> minB, maxB;
	vector<point<T> > valor;
	node<T>* hijo[8];
	node(point<T> a, point<T> b) : minB(a), maxB(b){
		for (int i = 0; i < 8; ++i) hijo[i] = nullptr;
	}
};


template<typename T, size_t CAPACIDAD>
class Octree{
	node<T>* root;

	bool fit(point<T> point, node<T> *node){
		bool X = point.x < (node->maxB).x && point.x >(node->minB).x;
		bool Y = point.y < (node->maxB).y && point.y >(node->minB).y;
		bool Z = point.z < (node->maxB).z && point.z >(node->minB).z;
		return (X && Y && Z);
	}

	bool find(point<T> val, node<T>* &ptr, node<T>* &par, node<T>* inicio = nullptr){
		node<T> *next;
		par = nullptr;
		ptr = inicio ? inicio : root;
		for (; ptr; par = ptr, ptr = next) {
			next = nullptr;
			for (int i = 0; i < 8; ++i){
				if (ptr->hijo[i] && fit(val, ptr->hijo[i])) next = ptr->hijo[i];
			}
			if (!next){
				for (auto point : ptr->valor)
					if (point == val) return true;
				return false;
			}
		}
		return false; 
	}

public:

	Octree(point<T> lb_left, point<T> uf_right){
		root = new node<T>(lb_left, uf_right);
	}

	bool find(point<T> point){
		node<T>* tmp, *tmp2;
		return (fit(point, root) && find(point, tmp, tmp2));
	}

	void insert(point<T> punto, node<T> *head = nullptr){
		node<T> *pos, *tmp;
		if (fit(punto, root) && !find(punto, pos, tmp, head)){
			if (pos->valor.size() == CAPACIDAD){
				vector<point<T> > tmp = pos->valor;
				tmp.push_back(punto);
				pos->valor.clear();
				T _x = pos->minB.x, _y = pos->minB.y, _z = pos->minB.z;
				T _X = pos->maxB.x, _Y = pos->maxB.y, _Z = pos->maxB.z;
				T _mx = (_x + _X) / 2, _my = (_y + _Y) / 2, _mz = (_z + _Z) / 2;
				pos->hijo[0] = new node<T>(point<T>(_x, _y, _z), point<T>(_mx, _my, _mz));
				pos->hijo[1] = new node<T>(point<T>(_mx, _y, _z), point<T>(_X, _my, _mz));
				pos->hijo[2] = new node<T>(point<T>(_x, _my, _z), point<T>(_mx, _Y, _mz));
				pos->hijo[3] = new node<T>(point<T>(_mx, _my, _z), point<T>(_X, _Y, _mz));
				pos->hijo[4] = new node<T>(point<T>(_x, _y, _mz), point<T>(_mx, _my, _Z));
				pos->hijo[5] = new node<T>(point<T>(_mx, _y, _mz), point<T>(_X, _my, _Z));
				pos->hijo[6] = new node<T>(point<T>(_x, _my, _mz), point<T>(_mx, _Y, _Z));
				pos->hijo[7] = new node<T>(point<T>(_mx, _my, _mz), point<T>(_X, _Y, _Z));
				for (auto p : tmp) insert(p, pos);
			}
			else
			{
				pos->valor.push_back(punto);
			}
		}
	}

	void erase(point<T> point){
		node<T> *pos, *parent;
		if (fit(point, root) && find(point, pos, parent)){
			auto it = pos->valor.begin();
			while (*it != point) ++it;
			pos->valor.erase(it);
			if (parent){
				int count = 0;
				for (int i = 0; i < 8; ++i) count += parent->hijo[i]->valor.size();
				if (count == CAPACIDAD){
					for (int i = 0; i < 8; ++i){
						for (auto p : parent->hijo[i]->valor)
							parent->valor.push_back(p);
						delete parent->hijo[i];
						parent->hijo[i] = nullptr;
					}
				}
			}
		}
	}

};

const float inf = 1e9;
enum { INSERT, ERASE };

struct tools {

	template<typename T>
	static pair<point<T>, point<T> > limites(char path[]){
		ifstream points(path);
		string line, tag;
		T x, y, z, mx, my, mz, Mx, My, Mz;
		mx = my = mz = inf;
		Mx = My = Mz = -inf;
		while (getline(points, line)){
			size_t pos = line.find(' ');
			if (pos == string::npos) continue;
			tag = line.substr(0, pos);
			if (tag != "v") continue;
			stringstream point(line.substr(pos));
			point >> x >> y >> z;
			mx = min(mx, x);
			my = min(my, y);
			mz = min(mz, z);
			Mx = max(Mx, x);
			My = max(My, y);
			Mz = max(Mz, z);
		}
		mx -= (Mx - mx) / 10;
		my -= (My - my) / 10;
		mz -= (Mz - mz) / 10;
		Mx += (Mx - mx) / 10;
		My += (My - my) / 10;
		Mz += (Mz - mz) / 10;
		cout << "Limite inferior : (" << mx << "," << my << "," << mz << ")\n";
		cout << "Limite superior : (" << Mx << "," << My << "," << Mz << ")\n";
		points.close();
		return make_pair(point<T>(mx, my, mz), point<T>(Mx, My, Mz));
	}

	template<typename T, int S>
	static int administrarP(Octree<T, S> &tree, int cod, char path[]){
		ifstream points(path);
		string line, tag;
		T x, y, z;
		int count = 0;
		while (getline(points, line)){
			size_t pos = line.find(' ');
			if (pos == string::npos) continue;
			tag = line.substr(0, pos);
			if (tag != "v") continue;
			stringstream punto(line.substr(pos));
			punto >> x >> y >> z;
			if (cod == INSERT) tree.insert(point<T>(x, y, z));
			if (cod == ERASE) tree.erase(point<T>(x, y, z));
			count++;
		}
		points.close();
		return count;
	}

};

const size_t capacidad = 5;
char archivo_obj[] = "../../archivo/cat.obj";

int main(int argc, char* args[]){
	char* objeto = (argc == 2 ? args[1] : archivo_obj);

	pair<point<float>, point<float> > bounds = tools::limites<float>(objeto);
	Octree<float, capacidad> tree(bounds.first, bounds.second);

	int cntI = tools::administrarP<float, capacidad>(tree, INSERT, objeto);
	cout << cntI << " inserted points" << endl;

	int cntO = tools::administrarP<float, capacidad>(tree, ERASE, objeto);
	cout << cntO << " deleted points" << endl;
}
