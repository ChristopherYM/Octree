#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <sstream>

#include <QtCore/QObject>
#include <Qt3DExtras/QCuboidMesh>
#include <Qt3DExtras/QPhongMaterial>
#include <QGuiApplication>
#include <Qt3DRender/qcamera.h>
#include <Qt3DCore/qentity.h>
#include <Qt3DRender/qcameralens.h>
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QCommandLinkButton>
#include <QtGui/QScreen>
#include <Qt3DExtras/qtorusmesh.h>
#include <Qt3DRender/qmesh.h>
#include <Qt3DRender/qtechnique.h>
#include <Qt3DRender/qmaterial.h>
#include <Qt3DRender/qeffect.h>
#include <Qt3DRender/qtexture.h>
#include <Qt3DRender/qrenderpass.h>
#include <Qt3DRender/qsceneloader.h>
#include <Qt3DRender/qpointlight.h>
#include <Qt3DCore/qtransform.h>
#include <Qt3DCore/qaspectengine.h>
#include <Qt3DRender/qrenderaspect.h>
#include <Qt3DExtras/qforwardrenderer.h>
#include <Qt3DExtras/qt3dwindow.h>
#include <Qt3DExtras/qfirstpersoncameracontroller.h>

using namespace std;

class SceneModifier : public QObject{
	Q_OBJECT

public:
	explicit SceneModifier(Qt3DCore::QEntity *rootEntity, float a, float b, float c, float d) {
		cuboid = new Qt3DExtras::QCuboidMesh();

		cuboidTransform = new Qt3DCore::QTransform();
		cuboidTransform->setScale(s);
		cuboidTransform->setTranslation(QVector3D(x, y, z));

		cuboidMaterial = new Qt3DExtras::QPhongMaterial();
		cuboidMaterial->setDiffuse(QColor(QRgb(0x665423)));

		m_cuboidEntity = new Qt3DCore::QEntity(m_rootEntity);
		m_cuboidEntity->addComponent(cuboid);
		m_cuboidEntity->addComponent(cuboidMaterial);
		m_cuboidEntity->addComponent(cuboidTransform);

	}
	~SceneModifier() {
		delete cuboid;
		delete cuboidTransform;
		delete cuboidMaterial;
	}

private:
	float x, y, z, s;
	Qt3DCore::QEntity *m_rootEntity;
	Qt3DCore::QEntity *m_cuboidEntity;
	Qt3DCore::QTransform *cuboidTransform;

	Qt3DExtras::QCuboidMesh *cuboid;
	Qt3DExtras::QPhongMaterial *cuboidMaterial;
};


enum { WHITE, BLACK, INTERNAL };

template<typename T>
struct point{
	T x, y, z;
	point(T a, T b, T c) : x(a), y(b), z(c) {}
};

template<typename T>
struct node{
	point<T> minB, maxB;
	int color;
	node<T>* hijo[8];
	node(point<T> a, point<T> b, int c) : minB(a), maxB(b), color(c){
		for (int i = 0; i < 8; ++i) hijo[i] = nullptr;
	}
};

template<typename T, size_t DEPTH>
class Octree{
	node<T>* root;

	bool fit(point<T> point, node<T> *node){
		bool X = point.x < (node->maxB).x && point.x >= (node->minB).x;
		bool Y = point.y < (node->maxB).y && point.y >= (node->minB).y;
		bool Z = point.z < (node->maxB).z && point.z >= (node->minB).z;
		return (X && Y && Z);
	}

	bool find(point<T> val, stack<node<T>* > &ptr, int col, node<T>* inicio = nullptr){
		node<T> *current, *next;
		for (current = inicio ? inicio : root; current; current = next) {
			next = nullptr;
			for (int i = 0; i < 8; ++i)
			{
				if (current->hijo[i] && fit(val, current->hijo[i])) next = current->hijo[i];
			}
			ptr.push(current);
		}
		return ptr.top()->color == col;
	}

public:

	Octree(point<T> lb_left, point<T> uf_right){
		root = new node<T>(lb_left, uf_right, WHITE);
	}

	bool find(point<T> point){
		stack<node<T>*> tmp;
		return (fit(point, root) && find(point, tmp, BLACK));
	}

	void insert(point<T> punto, int col = BLACK){
		stack<node<T>*> path;
		if (fit(punto, root) && !find(point, path, col)){
			while (path.size() != DEPTH){
				node<T>* pos = path.top();
				path.pop();
				int cc = pos->color;
				pos->color = INTERNAL;
				T _x = pos->minB.x, _y = pos->minB.y, _z = pos->minB.z;
				T _X = pos->maxB.x, _Y = pos->maxB.y, _Z = pos->maxB.z;
				T _mx = (_x + _X) / 2, _my = (_y + _Y) / 2, _mz = (_z + _Z) / 2;
				pos->hijo[0] = new node<T>(punto<T>(_x, _y, _z), point<T>(_mx, _my, _mz), cc);
				pos->hijo[1] = new node<T>(point<T>(_mx, _y, _z), point<T>(_X, _my, _mz), cc);
				pos->hijo[2] = new node<T>(point<T>(_x, _my, _z), point<T>(_mx, _Y, _mz), cc);
				pos->hijo[3] = new node<T>(point<T>(_mx, _my, _z), point<T>(_X, _Y, _mz), cc);
				pos->hijo[4] = new node<T>(point<T>(_x, _y, _mz), point<T>(_mx, _my, _Z), cc);
				pos->hijo[5] = new node<T>(point<T>(_mx, _y, _mz), point<T>(_X, _my, _Z), cc);
				pos->hijo[6] = new node<T>(point<T>(_x, _my, _mz), point<T>(_mx, _Y, _Z), cc);
				pos->hijo[7] = new node<T>(point<T>(_mx, _my, _mz), point<T>(_X, _Y, _Z), cc);
				find(punto, path, col, pos);
			}
			path.top()->color = col;
			path.pop();
			bool flag = true;
			while (!path.empty() && flag){
				node<T>* current = path.top();
				path.pop();
				for (int i = 0; i < 8; ++i){
					if (current->hijo[i]->color != col) flag = false;
				}
				if (flag){
					current->color = col;
					for (int i = 0; i < 8; ++i){
						delete current->hijo[i];
						current->hijo[i] = nullptr;
					}
				}
			}
		}
	}

	void erase(point<T> point){
		insert(point, WHITE);
	}

	void dfs(Qt3DCore::QEntity *entity, node<T> *cur = nullptr){
		if (cur == nullptr) cur = root;
		if (cur->hijo[0]){
			for (int i = 0; i < 8; ++i) dfs(entity, cur->hijo[i]);
		}
		else if (cur->color == BLACK){
			T cx, cy, cz, cube_size;
			cx = (cur->minB.x + cur->maxB.x) / 2;
			cy = (cur->minB.y + cur->maxB.y) / 2;
			cz = (cur->minB.z + cur->maxB.z) / 2;
			cube_size = cur->maxB.x - cur->minB.x;
			SceneModifier *p = new SceneModifier(entity, cx, cy, cz, cube_size);
		}
	}

};

const float inf = 1e9;
enum { INSERT, ERASE };

struct tools {
	template<typename T>
	static pair<point<T>, point<T> > GetBounds(char path[]){
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
		// Expand the bound by 10% and make a cube
		T dist = max(Mx - mx, max(My - my, Mz - mz));
		T plus = dist / 10;
		Mx = mx + dist + plus;
		My = my + dist + plus;
		Mz = mz + dist + plus;
		mx -= plus;
		my -= plus;
		mz -= plus;
		cout << "Limite inferior: (" << mx << "," << my << "," << mz << ")\n";
		cout << "Limite superior: (" << Mx << "," << My << "," << Mz << ")\n";
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

const size_t depth = 7;
char sample_path[] = "../../archivo/cat.obj";

void build_tree(Qt3DCore::QEntity *entity){
	pair<point<float>, point<float> > bounds = tools::GetBounds<float>(sample_path);
	Octree<float, depth> tree(bounds.first, bounds.second);

	int cntIn = tools::administrarP<float, depth>(tree, INSERT, sample_path);
	cout << cntIn << " inserted points" << endl;

	tree.dfs(entity);
}

int main(int argc, char **argv){
	QApplication app(argc, argv);
	Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
	Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity();

	view->defaultFrameGraph()->setClearColor(QColor(QRgb(0x4d4d4f)));

	QWidget *container = QWidget::createWindowContainer(view);
	QSize screenSize = view->screen()->size();

	QWidget *widget = new QWidget;
	QHBoxLayout *hLayout = new QHBoxLayout(widget);
	hLayout->addWidget(container, 1);

	widget->setWindowTitle(QStringLiteral("Octree visualization"));

	Qt3DRender::QCamera *cameraEntity = view->camera();
	cameraEntity->setPosition(QVector3D(35, 45, 30)); // Max Bound
	cameraEntity->setViewCenter(QVector3D(0, -10, -35)); // minBound

	Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(rootEntity);
	Qt3DRender::QPointLight *light = new Qt3DRender::QPointLight(lightEntity);
	light->setColor("white");
	light->setIntensity(3);
	lightEntity->addComponent(light);
	Qt3DCore::QTransform *lightTransform = new Qt3DCore::QTransform(lightEntity);
	lightTransform->setTranslation(cameraEntity->position());
	lightEntity->addComponent(lightTransform);

	Qt3DExtras::QFirstPersonCameraController *camController = new Qt3DExtras::QFirstPersonCameraController(rootEntity);
	camController->setCamera(cameraEntity);

	build_tree(rootEntity);

	view->setRootEntity(rootEntity);

	widget->show();
	widget->resize(1200, 800);

	return app.exec();
}
