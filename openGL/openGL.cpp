// openGL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
#include <iostream>

#include <GL/glut.h>
//#include <GL/glew.h>
//#include "Eigen/Dense"

using namespace std;
//using namespace Eigen;

//#pragma comment(lib, "glew32.lib")

const double TWO_PI = 6.2831853;
const double PI = 3.1415926;

GLsizei winWidth = 600, winHeight = 500;
int xRaster = 25, yRaster = 100;
GLubyte label[36] = { 'J', 'a', 'n',		'F', 'e', 'b',		'M', 'a', 'r',
					  'A', 'p', 'r',		'M', 'a', 'y',		'J', 'u', 'n',
					  'J', 'u', 'l',		'A', 'u', 'g',		'S', 'e', 'p',
					  'O', 'c', 't',		'N', 'o', 'v',		'D', 'e', 'c',
};

int dataValue[12] = { 420, 342, 324, 310, 262, 185, 190, 196, 217, 240, 312, 438 };

float xwcMin = 0.0, xwcMax = 225.0;
float ywcMin = 0.0, ywcMax = 225.0;

class screePt
{
private:
	int x;
	int y;

public:
	screePt()
	{
		x = y = 0;
	}

	void setCoords(int xValue, int yValue)
	{
		x = xValue;
		y = yValue;
	}

	int getX() const
	{
		return x;
	}

	int getY()const
	{
		return y;
	}

	void incrementX()
	{
		x++;
	}

	void decrementY()
	{
		y--;
	}
};

class wcPt2D
{	
public:
	float x, y;
	wcPt2D() { x = y = 0.0; }
	void setCoords(float vx, float vy) { x = vx; y = vy; }
	~wcPt2D() { /*cout << "destroy wcpt2D'x:" << x << " y:" << y << endl;*/ }
	//wcPt2D(wcPt2D&& value) { x = value.x; y = value.y; }
};

typedef float Matrix3x3[3][3];
Matrix3x3 matComposite;

class wcPt3D
{
public:
	float x, y, z, w;
	wcPt3D() { x = y = z = w = 0; }
	void setCoords(float xCoord, float yCoord, float zCoord, float wCoord = 1)
	{
		x = xCoord;
		y = yCoord;
		z = zCoord;
		w = wCoord;
	}

	void homogeneous()
	{
		x = x / w;
		y = y / w;
		w = 1;
	}
};

typedef float Matrix4x4[4][4];

inline int roundPlus(const float a) { return int(a + 0.5f); }
inline float dotPlus(const wcPt3D& p1, const wcPt3D& p2){ return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z; }
inline wcPt3D crossPlus(const wcPt3D& p1, const wcPt3D& p2)
{
	wcPt3D pt;
	pt.setCoords(p1.y*p2.z - p1.z*p2.y, p1.z*p2.x - p1.x*p2.z, p1.x*p2.y - p1.y*p2.x);
	return pt;
}
inline wcPt3D dir(const wcPt3D& p1, const wcPt3D& p2)
{
	wcPt3D result;
	result.setCoords(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	return result;
}
inline void normalize(wcPt3D& pt)
{
	float len = sqrt(dotPlus(pt, pt));
	pt.x/=len;
	pt.y/=len;
	pt.z/=len;
}

void setPixel(GLint x, GLint y)
{
	glBegin(GL_POINTS);
	glVertex2i(x, y);
	glEnd();
}

void swapPts(wcPt2D* p1, wcPt2D* p2)
{
	wcPt2D temp;
	temp = *p1;
	*p1 = *p2;
	*p2 = temp;
}

#pragma region Matrix3x3
//
void matrix3x3SetIdentity(Matrix3x3& matIdent3x3)
{
	int row, col;
	for (row = 0; row < 3; ++row)
	{
		for (col = 0; col < 3; ++col)
		{
			matIdent3x3[row][col] = (row == col);
		}
	}
}

void matrix3x3PreMultiply(Matrix3x3& m1, Matrix3x3& m2)
{
	int row, col;
	Matrix3x3 matTemp;
	for (row = 0; row<3;++row)
	{
		for (col = 0; col < 3; ++col)
		{
			matTemp[row][col] = m1[row][0] * m2[0][col] + m1[row][1] * m2[1][col] + m1[row][2] * m2[2][col];
		}
	}

	for (row = 0; row < 3; ++row)
	{
		for (col = 0; col < 3; ++col)
		{
			m2[row][col] = matTemp[row][col];
		}
	}
}

void translate2D(float tx, float ty)
{
	Matrix3x3 matTrans;
	matrix3x3SetIdentity(matTrans);

	matTrans[0][2] = tx;
	matTrans[1][2] = ty;

	matrix3x3PreMultiply(matTrans, matComposite);
}

void rotate2D(wcPt2D& pivotPt, float theta)
{
	Matrix3x3 matRot;
	matrix3x3SetIdentity(matRot);

	matRot[0][0] = cos(theta);
	matRot[0][1] = -sin(theta);
	matRot[0][2] = pivotPt.x * (1 - cos(theta)) + pivotPt.y * sin(theta);

	matRot[1][0] = sin(theta);
	matRot[1][1] = cos(theta);
	matRot[1][2] = pivotPt.y * (1 - cos(theta)) - pivotPt.x * sin(theta);

	matrix3x3PreMultiply(matRot, matComposite);
}

void scale2D(float sx, float sy, wcPt2D& fixedPt)
{
	Matrix3x3 matScale;
	matrix3x3SetIdentity(matScale);

	matScale[0][0] = sx;
	matScale[0][2] = (1 - sx) * fixedPt.x;

	matScale[1][1] = sy;
	matScale[1][2] = (1 - sy) * fixedPt.y;

	matrix3x3PreMultiply(matScale, matComposite);
}

void transformVerts2D(int nVerts, wcPt2D* verts)
{
	int k;
	float temp;

	for (k = 0; k < nVerts; ++k)
	{
		temp = matComposite[0][0] * verts[k].x + matComposite[0][1] * verts[k].y + matComposite[0][2];
		verts[k].y = matComposite[1][0] * verts[k].x + matComposite[1][1] * verts[k].y + matComposite[1][2];
		verts[k].x = temp;
	}
}

#pragma endregion

#pragma region Matrix4x4
Matrix4x4 matComposite3D;
void matrix4x4SetIdentity(Matrix4x4& matIdentity)
{
	int row, col;
	for (row = 0; row < 4; ++row)
		for (col = 0; col < 4; ++col)
			matIdentity[row][col] = col == row;
}

void matrix4x4PreMultiply(Matrix4x4& m1, Matrix4x4& m2)
{
	int row, col;
	Matrix4x4 matTemp;
	for (row = 0; row < 4; ++row)
	{
		for (col = 0; col < 4; ++col)
		{
			matTemp[row][col] = m1[row][0] * m2[0][col] + m1[row][1] * m2[1][col] + m1[row][2] * m2[2][col] + m1[row][3]*m2[3][col];
		}
	}

	for (row = 0; row < 4; ++row)
	{
		for (col = 0; col < 4; ++col)
		{
			m2[row][col] = matTemp[row][col];
		}
	}
}

//
void translate3D(float tx, float ty, float tz)
{
	Matrix4x4 matTrans3D;
	matrix4x4SetIdentity(matTrans3D);

	matTrans3D[0][3] = tx;
	matTrans3D[1][3] = ty;
	matTrans3D[2][3] = tz;

	matrix4x4PreMultiply(matTrans3D, matComposite3D);
}

//
void rotate3D(const wcPt3D& p1, const wcPt3D& p2, float radianAngle)
{
	Matrix4x4 matQuaternionRot;
	float axisVectLen = sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z));

	float sinA = sin(radianAngle);
	float cosA = cos(radianAngle);
	float oneC = 1 - cosA;

	float ux = (p2.x - p1.x) / axisVectLen;
	float uy = (p2.y - p1.y) / axisVectLen;
	float uz = (p2.z - p1.z) / axisVectLen;

	translate3D(-p1.x, -p1.y, -p1.z);
	matrix4x4SetIdentity(matQuaternionRot);
	matQuaternionRot[0][0] = ux * ux * oneC + cosA;
	matQuaternionRot[0][1] = ux * uy * oneC - uz * sinA;
	matQuaternionRot[0][2] = ux * uz * oneC + uy * sinA;
	matQuaternionRot[1][0] = uy * ux * oneC + uz * sinA;
	matQuaternionRot[1][1] = uy * uy * oneC + cosA;
	matQuaternionRot[1][2] = uy * uz * oneC - ux * sinA;
	matQuaternionRot[2][0] = uz * ux * oneC - uy * sinA;
	matQuaternionRot[2][1] = uz * uy * oneC + ux * sinA;
	matQuaternionRot[2][2] = uz * uz * oneC + cosA;

	matrix4x4PreMultiply(matQuaternionRot, matComposite3D);

	translate3D(p1.x, p1.y, p1.z);
}

void scale3D(float sx, float sy, float sz, const wcPt3D& fixedPt)
{
	Matrix4x4 matScale3D;
	matrix4x4SetIdentity(matScale3D);

	matScale3D[0][0] = sx;
	matScale3D[0][3] = (1 - sx) * fixedPt.x;
	matScale3D[1][0] = sy;
	matScale3D[1][3] = (1 - sy) * fixedPt.y;
	matScale3D[2][0] = sz;
	matScale3D[2][3] = (1 - sz) * fixedPt.z;

	matrix4x4PreMultiply(matScale3D, matComposite3D);
}

void transformVerts3D(int nVerts, wcPt3D* verts)
{
	int k;
	for (k = 0; k < nVerts; ++k)
	{
		float tempX = verts[k].x * matComposite3D[0][0] + verts[k].y * matComposite3D[0][1] + verts[k].z * matComposite3D[0][2] + verts[k].w * matComposite3D[0][3];
		float tempY = verts[k].x * matComposite3D[1][0] + verts[k].y * matComposite3D[1][1] + verts[k].z * matComposite3D[1][2] + verts[k].w * matComposite3D[1][3];
		float tempZ = verts[k].x * matComposite3D[2][0] + verts[k].y * matComposite3D[2][1] + verts[k].z * matComposite3D[2][2] + verts[k].w * matComposite3D[2][3];
		verts[k].w = verts[k].x * matComposite3D[3][0] + verts[k].y * matComposite3D[3][1] + verts[k].z * matComposite3D[3][2] + verts[k].w * matComposite3D[3][3];
		verts[k].z = tempZ;
		verts[k].y = tempY;
		verts[k].x = tempX;
	}
}

#pragma endregion

#pragma region OpenGL Draw Sampler

void lineSegment(screePt p1, screePt p2)
{
	glBegin(GL_LINES);
	glVertex2i(p1.getX(), p1.getY());
	glVertex2i(p2.getX(), p2.getY());
	glEnd();
}

void triangle2D(wcPt2D* verts, GLint nVerts)
{
	GLint k;

	glPolygonMode(GL_FRONT, GL_FILL);
	glBegin(GL_POLYGON);
	for (k = 0; k < nVerts; ++k)
	{
		glVertex2f(verts[k].x, verts[k].y);
	}
	glEnd();
}

void triangle3D(wcPt3D* verts, GLint nVerts)
{
	GLint k;

	glPolygonMode(GL_FRONT, GL_FILL);
	glBegin(GL_POLYGON);
	for (k = 0; k < nVerts; ++k)
	{
		glVertex3f(verts[k].x, verts[k].y, 0);//此处先不管Z值
	}
	glEnd();
}

typedef GLint vertex3[3];
vertex3 pt[8] = { {0, 0, 0}, {0, 100, 0}, {100, 0, 0}, {100, 100, 0}, {0, 0, 100}, {0, 100, 100}, {100, 0, 100}, {100, 100, 100} };
void quad(GLint n1, GLint n2, GLint n3, GLint n4)
{
	glBegin(GL_QUADS);
	glVertex3iv(pt[n1]);
	glVertex3iv(pt[n2]);
	glVertex3iv(pt[n3]);
	glVertex3iv(pt[n4]);
	glEnd();
}

void cube()
{
	quad(6, 2, 3, 7);
	quad(5, 1, 0, 4);
	quad(7, 3, 1, 5);
	quad(4, 0, 2, 6);
	quad(2, 0, 1, 3);
	quad(4, 6, 7, 5);
}

void cubeElements()
{
	glEnableClientState(GL_VERTEX_ARRAY);//激活顶点数组
	glVertexPointer(3, GL_INT, 0, pt);//设置顶点坐标的位置和格式

	GLubyte vertexIndex[] = { 6, 2, 3, 7, 5, 1, 0, 4, 7, 3, 1, 5, 4, 0, 2, 6, 2, 0, 1, 3, 7, 5, 4, 6 };//设置顶点索引
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_BYTE, vertexIndex);
}

void callList()
{
	GLuint regHex;
	GLdouble theta;
	GLint x, y, k;

	regHex = glGenLists(1);
	glNewList(regHex, GL_COMPILE);
	glBegin(GL_POLYGON);
	for (k = 0; k < 6; ++k)
	{
		theta = TWO_PI * k / 6.0;
		x = 200 + 150 * cos(theta);
		y = 200 + 150 * sin(theta);
		glVertex2i(x, y);
	}
	glEnd();
	glEndList();

	glCallList(regHex);
}

void lineGraph()
{
	GLint month, k;
	GLint x = 30;

	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_LINE_STRIP);
	for (k = 0; k < 12; k++)
	{
		glVertex2i(x + k * 50, dataValue[k]);
	}
	glEnd();

	glColor3f(1.0f, 0.0f, 0.0f);
	xRaster = 25;
	for (k = 0; k < 12; k++)
	{
		glRasterPos2i(xRaster + k * 50, dataValue[k] - 4);
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '*');
	}

	glColor3f(0.0f, 0.0f, 0.0f);
	xRaster = 20;
	for (month = 0; month < 12; month++)
	{
		glRasterPos2i(xRaster, yRaster);
		for (k = 3 * month; k < 3 * month + 3; k++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[k]);
		}
		xRaster += 50;
	}
}

void barChart()
{
	GLint month, k;

	glColor3f(1.0, 0.0, 0.0);
	for (k = 0; k < 12; ++k)
	{
		glRecti(20 + 50 * k, 150, 40 + 50 * k, dataValue[k]);
	}

	glColor3f(0.0f, 0.0f, 0.0f);
	xRaster = 20;
	for (month = 0; month < 12; month++)
	{
		glRasterPos2i(xRaster, yRaster);
		for (k = 3 * month; k < 3 * month + 3; k++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[k]);
		}
		xRaster += 50;
	}
}

typedef enum
{
	limacon = 1,
	cardioid,
	threeLeaf,
	fourLeaf,
	spiral,
} curveNames;

GLuint drawCurve(int curveNum)
{
	GLuint regHex = glGenLists(1);
	glNewList(regHex, GL_COMPILE);

	const int a = 175, b = 60;

	float r, theta, dtheta = 1.0 / (float)(a);
	int x0 = 200, y0 = 250;

	screePt curvePt[2];

	glColor3f(0.0f, 0.0f, 0.0f);

	curvePt[0].setCoords(x0, y0);

	switch (curveNum)
	{
	case limacon:
		curvePt[0].setCoords(curvePt[0].getX() + a + b, curvePt[0].getY());
		break;
	case cardioid:
		curvePt[0].setCoords(curvePt[0].getX() + a + a, curvePt[0].getY());
		break;
	case threeLeaf:
		curvePt[0].setCoords(curvePt[0].getX() + a, curvePt[0].getY());
		break;
	case fourLeaf:
		curvePt[0].setCoords(curvePt[0].getX() + a, curvePt[0].getY());
		break;
	case spiral:
		break;
	default:
		break;
	}

	theta = dtheta;
	while (theta < TWO_PI)
	{
		switch (curveNum)
		{
		case limacon:
			r = a * cos(theta) + b; break;
		case cardioid:
			r = a * (1 + cos(theta)); break;
		case threeLeaf:
			r = a * cos(3 * theta); break;
		case fourLeaf:
			r = a * cos(2 * theta); break;
		case spiral:
			r = (a / 4.0) * theta; break;
		default:
			break;
		}

		curvePt[1].setCoords(x0 + r * cos(theta), y0 + r * sin(theta));
		lineSegment(curvePt[0], curvePt[1]);

		curvePt[0].setCoords(curvePt[1].getX(), curvePt[1].getY());
		theta += dtheta;
	}

	glEndList();
	return regHex;
}

void initDrawList()
{
	for (int i = 1; i < 6; ++i)
	{
		drawCurve(i);
	}
}

void drawColorTriangles()
{
	glShadeModel(GL_SMOOTH);//GL_SMOOTH-平滑插值，GL_FLAT-单色填充
	glPolygonMode(GL_FRONT, GL_FILL);//GL_LINE-显示边框，GL_POINT-显示顶点，GL_FILL-默认填充显示
	glEnable(GL_LINE_SMOOTH);

	glBegin(GL_TRIANGLES);
	glColor3f(.0f, .0f, 1.f);
	glVertex2i(50, 50);
	glColor3f(.0f, 1.f, .0f);
	glVertex2i(250, 50);
	glColor3f(1.f, .0f, .0f);
	glVertex2i(150, 250);
	glEnd();
}

#pragma endregion

#pragma region Raster
void lineDDA(int x0, int y0, int xEnd, int yEnd)
{
	float dx = xEnd - x0;
	float dy = yEnd - y0;
	int steps, k;

	float xIncrement, yIncrement, x = x0, y = y0;

	if (abs(dx) > abs(dy))
	{
		steps = dx;
	}
	else
	{
		steps = dy;
	}

	xIncrement = dx / steps;
	yIncrement = dy / steps;

	setPixel(roundPlus(x), roundPlus(y));
	for (k = 0; k < steps; ++k)
	{
		x = x0 + xIncrement * k;
		y = y0 + yIncrement * k;
		setPixel(roundPlus(x), roundPlus(y));
	}
}

void lineBres(int x0, int y0, int xEnd, int yEnd)
{
	if (xEnd < x0)
		swap(xEnd, x0);
	if (yEnd < y0)
		swap(yEnd, y0);

	int dx = xEnd - x0;
	int dy = yEnd - y0;
	int dx2 = dx + dx;
	int dy2 = dy + dy;
	int x, y;

	if (x0 > xEnd)
	{
		x = xEnd;
		y = yEnd;
	}
	else
	{
		x = x0;
		y = y0;
	}

	setPixel(x, y);
	if (dx > dy)
	{
		int p = dy2 - dx;
		int a = dy2 - dx2;
		for (int k = 0; k < dx; ++k)
		{
			x += 1;
			if (p < 0)
			{
				p = p + dy2;
			}
			else
			{
				y += 1;
				p = p + a;
			}
			setPixel(x, y);
		}
	} 
	else
	{
		int p = dx2 - dy;
		int a = dx2 - dy2;
		for (int k = 0; k < dy; ++k)
		{
			y += 1;
			if (p < 0)
			{
				p = p + dx2;
			}
			else
			{
				x += 1;
				p = p + a;
			}
			setPixel(x, y);
		}
	}
}

//中点法画圆
void circlePlotPoints(int x, int y, screePt point)
{
	setPixel(x + point.getX(), y + point.getY());
	setPixel(x + point.getX(), y - point.getY());
	setPixel(x - point.getX(), y + point.getY());
	setPixel(x - point.getX(), y - point.getY());
	setPixel(x + point.getY(), y + point.getX());
	setPixel(x + point.getY(), y - point.getX());
	setPixel(x - point.getY(), y + point.getX());
	setPixel(x - point.getY(), y - point.getX());
}
void circleMidPoint(int x, int y, GLint radius)
{
	screePt circPt;
	int p = 1 - radius;
	circPt.setCoords(0, radius);
	circlePlotPoints(x, y, circPt);

	while (circPt.getX() < circPt.getY())
	{
		circPt.incrementX();

		if (p < 0)
		{
			p += 2 * circPt.getX() + 1;
		}
		else
		{
			circPt.decrementY();
			p += 2 * circPt.getX() - 2 * circPt.getY() + 1;
		}

		circlePlotPoints(x, y, circPt);
	}
}

//中点画椭圆
void ellipsePlotPoints(int xCenter, int yCenter, int x, int y)
{
	setPixel(xCenter + x, yCenter + y);
	setPixel(xCenter + x, yCenter - y);
	setPixel(xCenter - x, yCenter + y);
	setPixel(xCenter - x, yCenter - y);
}
void ellipseMidPoint(int xCenter, int yCenter, int xRadius, int yRadius)
{
	const int Rx2 = xRadius * xRadius;
	const int Ry2 = yRadius * yRadius;
	const int tRx2 = 2 * Rx2;
	const int tRy2 = 2 * Ry2;

	int p = 0;
	int x = 0;
	int y = yRadius;
	int px = 0;
	int py = tRx2 * y;

	ellipsePlotPoints(xCenter, yCenter, x, y);

	p = roundPlus(Ry2 - (Rx2 * yRadius) + (0.25f * Rx2));
	while (px < py)
	{
		x++;
		px += tRy2;
		if (p < 0)
		{
			p += Ry2 + px;
		}
		else
		{
			y--;
			py -= tRx2;
			p += Ry2 + px - py;
		}
		ellipsePlotPoints(xCenter, yCenter, x, y);
	}

	p = roundPlus(Ry2 * (x + 0.5) * (x + 0.5) + Rx2 * (y - 1) * (y - 1) - Rx2 * Ry2);
	while (y > 0)
	{
		y--;
		py -= tRx2;
		if (p > 0)
		{
			p += Rx2 - py;
		}
		else
		{
			x++;
			px += tRy2;
			p += Rx2 - py + px;
		}
		ellipsePlotPoints(xCenter, yCenter, x, y);
	}
}

//中点抛物线
void midPointGravity(int vx0, int vy0, int x0, int y0)
{
	const int gravity = 980;
	const int vxy = vx0 * vy0;
	const int gx = gravity * (x0 - 0.5f);
	const int vxt = vx0 * vx0;
	const int xEnd = (vxy / gravity) + x0;

	int x = x0, y = y0;
	int xg = x * gravity;
	setPixel(x, y);
	setPixel(xEnd + xEnd - x, y);

	int temp = roundPlus((vy0 * vy0 - vxt) / (2 * gravity)) + y0;
	int p = vxt + 0.125 * gravity - 0.5 * vxy;
	while (y < temp)
	{
		y++;
		if (p > 0)
		{
			x++;
			xg += gravity;
			p += vxt + xg - gx - gravity - vxy;
		}
		else
		{
			p += vxt;
		}

		setPixel(x, y);
		setPixel(xEnd + xEnd - x, y);
	}

	p = (0.5f + y - y0) * vxt - vxy * (x + 1 - x0) + 0.5f * gravity * (x + 1 - x0) * (x + 1 - x0);
	while (x < xEnd)
	{
		x++;
		xg += gravity;
		if (p < 0)
		{
			y++;
			p += xg - gx - vxy + vxt;
		}
		else
		{
			p += xg - gx - vxy;
		}
		setPixel(x, y);
		setPixel(xEnd + xEnd - x, y);
	}
}

//绘制下平底三角形
void drawBottomFlatTrangle(int x0, int y0, int x1, int y1, int x2, int y2)
{
	if (y0 != y1)
	{
		cout << "this trangle is not bottom flat" << endl;
		return;
	}

	float deltax = x2 - x0;
	float deltay = y2 - y0;
	float deltaL = deltax / deltay;
	deltax = x2 - x1;
	deltay = y2 - y1;
	float deltaR = deltax / deltay;
	int k;

	for (k = y0; k < y2; ++k)
	{
		int xl = roundPlus((k - y0) * deltaL + x0);
		int xr = roundPlus((k - y1) * deltaR + x1);
		lineBres(xl, k, xr, k);
	}
}
//绘制上平底三角形
void drawTopFlatTrangle(int x0, int y0, int x1, int y1, int x2, int y2)
{
	if (y0 != y1)
	{
		cout << "this trangle is not top flat" << endl;
		return;
	}

	float deltax = x2 - x0;
	float deltay = y2 - y0;
	float deltaL = deltax / deltay;
	deltax = x2 - x1;
	deltay = y2 - y1;
	float deltaR = deltax / deltay;
	int k;

	for (k = y2; k < y0; ++k)
	{
		int xl = roundPlus((k - y0) * deltaL + x0);
		int xr = roundPlus((k - y1) * deltaR + x1);
		lineBres(xl, k, xr, k);
	}
}
//绘制任意三角形
void drawTrangle(int x0, int y0, int x1, int y1, int x2, int y2)
{
	//todo y2>y0>y1
	if (y0 < y1)
	{
		swap(y1, y0);
		swap(x1, x0);
	}
	if (y2 < y1)
	{
		swap(y2, y1);
		swap(x2, x1);
	}
	if (y2 < y0)
	{
		swap(y2, y0);
		swap(x2, x0);
	}

	if (y0 == y1)//下平底
		drawBottomFlatTrangle(x0, y0, x1, y1, x2, y2);
	else if (y0 == y2)//上平底
		drawTopFlatTrangle(x0, y0, x2, y2, x1, y1);
	else
	{
		float deltaX = x2 - x1;
		float dletaY = y2 - y1;
		float delta = deltaX / dletaY;
		int xAdd = roundPlus((y0 - y1) * delta + x1);
		int yAdd = y0;
		if (xAdd < x0)
		{
			swap(xAdd, x0);
			swap(yAdd, y0);
		}

		drawBottomFlatTrangle(x0, y0, xAdd, yAdd, x2, y2);
		drawTopFlatTrangle(x0, y0, xAdd, yAdd, x1, y1);
	}
}

void drawTrangle(const screePt& p1, const screePt& p2, const screePt& p3)
{
	drawTrangle(p1.getX(), p1.getY(), p2.getX(), p2.getY(), p3.getX(), p3.getY());
}

#pragma endregion

#pragma region 2D Clip

void drawClipWin(const wcPt2D& winMin, const wcPt2D& winMax)
{
	lineBres(winMin.x, winMin.y, winMax.x, winMin.y);
	lineBres(winMin.x, winMin.y, winMin.x, winMax.y);
	lineBres(winMax.x, winMin.y, winMax.x, winMax.y);
	lineBres(winMin.x, winMax.y, winMax.x, winMax.y);
}

const int winLeftBitCode = 0x1;
const int winRightBitCode = 0x2;
const int winBottomBitCode = 0x4;
const int winTopBitCode = 0x8;

inline int inside(int code) { return int(!code); }//取反
inline int reject(int code1, int code2) { return int(code1 & code2); }//且
inline int accept(int code1, int code2) { return int(!(code1 | code2)); }

GLubyte encode(wcPt2D& pt, wcPt2D& winMin, wcPt2D& winMax)
{
	GLubyte code = 0x00;
	if (pt.x < winMin.x)
		code = code | winLeftBitCode;
	if (pt.x > winMax.x)
		code |= winRightBitCode;
	if (pt.y < winMin.y)
		code |= winBottomBitCode;
	if (pt.y > winMax.y)
		code |= winTopBitCode;
	return code;
}

void swapCodes(GLubyte* c1, GLubyte* c2)
{
	GLubyte temp;
	temp = *c1;
	*c1 = *c2;
	*c2 = temp;
}

void lineClipCohSuth(wcPt2D& winMin, wcPt2D& winMax, wcPt2D& p1, wcPt2D& p2)
{
	GLubyte code1, code2;
	int done = false, plotLine = false;

	float m;

	while (!done)
	{
		code1 = encode(p1, winMin, winMax);
		code2 = encode(p2, winMin, winMax);

		if (accept(code1, code2))
		{
			done = true;
			plotLine = true;
		}
		else
		{
			if (reject(code1, code2))//同侧剔除
			{
				done = true;
			}
			else
			{
				if (inside(code1))
				{
					swapPts(&p1, &p2);
					swapCodes(&code1, &code2);
				}

				if (p2.x != p1.x)
					m = (p2.y - p1.y) / (p2.x - p1.x);
				if (code1 & winLeftBitCode)
				{
					p1.y += (winMin.x - p1.x) * m;
					p1.x = winMin.x;
				}
				else
				{
					if (code1 & winRightBitCode)
					{
						p1.y += (winMax.x - p1.x) * m;
						p1.x = winMax.x;
					}
					else
					{
						if (code1 & winBottomBitCode)
						{
							p1.x += (winMin.y - p1.y) / m;
							p1.y = winMin.y;
						}
						else
						{
							if (code1 & winTopBitCode)
							{
								p1.x += (winMax.y - p1.y) / m;
								p1.y = winMax.y;
							}
						}
					}
				}
			}
		}
	}

	if (plotLine)
		lineBres(roundPlus(p1.x), roundPlus(p1.y), roundPlus(p2.x), roundPlus(p2.y));
}

int clipTest(float p, float q, float* uIn, float* uOut)
{
	float r;
	int returnValue = true;

	if (p < .0)//入边取最大值
	{
		r = q / p;
		if (r > *uOut)
			returnValue = false;
		else
		{
			if (r > *uIn)
				*uIn = r;
		}
	}
	else
	{
		if (p > .0)//出边取最小值
		{
			r = q / p;
			if (r < *uIn)
				returnValue = false;
			else
			{
				if (r < *uOut)
					*uOut = r;
			}
		}
		else//与边界平行
			if (q < .0)
				returnValue = false;
	}
	return returnValue;
}

void lineClipLiangBarsk(wcPt2D& winMin, wcPt2D& winMax, wcPt2D& p1, wcPt2D& p2)
{
	float uIn = 0, uOut = 1, dx = p2.x - p1.x, dy;
	if (clipTest(-dx, p1.x - winMin.x, &uIn, &uOut))
	{
		if (clipTest(dx, winMax.x - p1.x, &uIn, &uOut))
		{
			dy = p2.y - p1.y;
			if (clipTest(-dy, p1.y - winMin.y, &uIn, &uOut))
			{
				if (clipTest(dy, winMax.y - p1.y, &uIn, &uOut))
				{
					if (uOut < 1.0)
						p2.setCoords(p1.x + uOut * dx, p1.y + uOut * dy);
					if (uIn > .0)
						p1.setCoords(p1.x + uIn * dx, p1.y + uIn * dy);

					lineBres(roundPlus(p1.x), roundPlus(p1.y), roundPlus(p2.x), roundPlus(p2.y));
				}
			}
		}
	}
}

void testLineClip()
{
	wcPt2D winMin, winMax;
	winMin.x = 0; winMin.y = 0;
	winMax.x = 100; winMax.y = 100;
	drawClipWin(winMin, winMax);

	wcPt2D p1, p2;
	p1.x = -50; p1.y = -50;
	p2.x = 120, p2.y = 120;
	lineClipCohSuth(winMin, winMax, p1, p2);
}

typedef enum Boundary
{
	Left,
	Right,
	Bottom,
	Top,
};
const int nClip = 4;

//是否在边界区域内
int inside(wcPt2D& p, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax)
{
	switch (edeg)
	{
	case Left:
		if (p.x < winMin.x)return (false);
		break;
	case Right:
		if (p.x > winMax.x)return (false);
		break;
	case Bottom:
		if (p.y < winMin.y)return (false);
		break;
	case Top:
		if (p.y > winMax.y)return (false);
		break;
	}
	return (true);
}

//是否穿过边界
int cross(wcPt2D& p1, wcPt2D& p2, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax)
{
	if (inside(p1, edeg, winMin, winMax) == inside(p2, edeg, winMin, winMax))
		return false;
	else 
		return true;
}

//获取边界交点
wcPt2D intersect(wcPt2D& p1, wcPt2D& p2, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax)
{
	wcPt2D iPt;
	float m;

	if (p1.x != p2.x) m = (p2.y - p1.y) / (p2.x - p1.x);
	switch (edeg)
	{
	case Left:
		iPt.x = winMin.x;
		iPt.y = p2.y + (winMin.x - p2.x) * m;
		break;
	case Right:
		iPt.x = winMax.x;
		iPt.y = p2.y + (winMax.x - p2.x) * m;
		break;
	case Bottom:
		iPt.y = winMin.y;
		if (p1.x != p2.x) iPt.x = p2.x + (winMin.y - p2.y) / m;
		else iPt.x = p2.x;
		break;
	case Top:
		iPt.y = winMax.y;
		if (p1.x != p2.x)iPt.x = p2.x + (winMax.y - p2.y) / m;
		else iPt.x = p2.x;
		break;
	default:
		break;
	}
	return iPt;
}

void clipPoint(wcPt2D& p, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax, wcPt2D* pOut, int* cnt, wcPt2D* first[], wcPt2D* last)
{
	wcPt2D iPt;
	if (!first[edeg])
	{
		first[edeg] = new wcPt2D();
		first[edeg]->setCoords(p.x, p.y);
	}
	else
	{
		if (cross(p, last[edeg], edeg, winMin, winMax))//与边界相交
		{
			iPt = intersect(p, last[edeg], edeg, winMin, winMax);
			if (edeg < Top)
				clipPoint(iPt, Boundary(edeg + 1), winMin, winMax, pOut, cnt, first, last);
			else
			{
				pOut[*cnt] = iPt;
				(*cnt)++;
			}
		}
	}

	last[edeg] = p;

	if (inside(p, edeg, winMin, winMax))//边界以内的点
	{
		if (edeg < Top)
			clipPoint(p, Boundary(edeg + 1), winMin, winMax, pOut, cnt, first, last);
		else
		{
			pOut[*cnt] = p;
			(*cnt)++;
		}
	}
}

void closeClip(wcPt2D& winMin, wcPt2D& winMax, wcPt2D* pOut, GLint* cnt, wcPt2D* first[], wcPt2D* last)
{
	wcPt2D pt;
	int edge;
	for (edge = Left; edge <= Top; ++edge)
	{
		if (cross(last[edge], *first[edge], Boundary(edge), winMin, winMax))
		{
			pt = intersect(last[edge], *first[edge], Boundary(edge), winMin, winMax);
			if (edge < Top)
				clipPoint(pt, Boundary(edge+1), winMin, winMax, pOut, cnt, first, last);
			else
			{
				pOut[*cnt] = pt;
				(*cnt)++;
			}
		}
	}
}

int polygonClipSuthHodg(wcPt2D& winMin, wcPt2D& winMax, GLint n, wcPt2D* pIn, wcPt2D* pOut)
{
	wcPt2D* first[nClip] = { 0, 0, 0, 0 };
	wcPt2D last[nClip];
	int k, cnt = 0;

	for (k = 0; k < n; ++k)
		clipPoint(pIn[k], Left, winMin, winMax, pOut, &cnt, first, last);

	closeClip(winMin, winMax, pOut, &cnt, first, last);
	return cnt;
}

void testPolygonClip()
{
	wcPt2D winMin, winMax;
	winMin.x = 0; winMin.y = 0;
	winMax.x = 100; winMax.y = 100; 
	drawClipWin(winMin, winMax);

	int nVerts = 3;
	wcPt2D p1, p2, p3;
	p1.setCoords(-25.0, 10.0);
	p2.setCoords(125.0, 10.0);
	p3.setCoords(50.0, 125.0);
	wcPt2D verts[3] = { p1, p2, p3};

	glColor3f(1.0f, 0.0f, 0.0f);
	triangle2D(verts, nVerts);

	wcPt2D outVerts[6];
	int outCnt = polygonClipSuthHodg(winMin, winMax, nVerts, verts, outVerts);
	glColor3f(0.0f, 1.0f, 0.0f);
	triangle2D(outVerts, outCnt);
}

#pragma endregion

#pragma region 3D View
//相机矩阵
void generateCameraModel(Matrix4x4& cameraMatrix, const wcPt3D& origin, const wcPt3D& lookAt, const wcPt3D& upDir)
{
	matrix4x4SetIdentity(cameraMatrix);

	wcPt3D zDir = dir(origin, lookAt);
	normalize(zDir);
	wcPt3D xDir = crossPlus(upDir, zDir);
	normalize(xDir);
	wcPt3D yDir = crossPlus(zDir, xDir);
	normalize(yDir);

	cameraMatrix[0][0] = xDir.x;
	cameraMatrix[0][1] = xDir.y;
	cameraMatrix[0][2] = xDir.z;
	cameraMatrix[0][3] = -dotPlus(origin, xDir);
	cameraMatrix[1][0] = yDir.x;
	cameraMatrix[1][1] = yDir.y;
	cameraMatrix[1][2] = yDir.z;
	cameraMatrix[1][3] = -dotPlus(origin, yDir);
	cameraMatrix[2][0] = zDir.x;
	cameraMatrix[2][1] = zDir.y;
	cameraMatrix[2][2] = zDir.z;
	cameraMatrix[2][3] = -dotPlus(origin, zDir);
}
//投影矩阵
void generateProjectModel(Matrix4x4& proMatrix, float fov, float aspect, float znear, float zfar)
{
	matrix4x4SetIdentity(proMatrix);

	float cot = atan(fov*0.5);

	proMatrix[0][0] = cot/aspect;
	proMatrix[1][1] = cot;
	proMatrix[2][2] = 0;//(znear + zfar)/(znear - zfar);
	proMatrix[2][3] = 0;//-(2*znear*zfar)/(znear - zfar);
	proMatrix[3][2] = 1;
	proMatrix[3][3] = 0;
}
//投屏矩阵
void generateScreenModel(Matrix4x4& screenMatrix, const wcPt3D& center, float width, float height)
{
	matrix4x4SetIdentity(screenMatrix);

	screenMatrix[0][0] = width / 2;
	screenMatrix[0][3] = center.x*width / 2;
	screenMatrix[1][1] = height / 2;
	screenMatrix[1][3] = center.y*height / 2;
	screenMatrix[2][2] = 0;
	screenMatrix[3][3] = 0;
}
#pragma endregion

//opengl api draw
void drawFunc(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.4f, 0.2f);//设定显示对象的颜色

	//lineSegment();
	//cube();
	//cubeElements();
	//callList();
	//lineGraph();
	//barChart();
	drawColorTriangles();

	glFlush();
}

// custom raster draw
void rasterFunc()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.4f, 0.2f);//设定显示对象的颜色

	//lineBres(20, 10, 330, 218);
	//circleMidPoint(200, 200, 100);
	//ellipseMidPoint(200, 200, 100, 50);
	//midPointGravity(200, 800, 40, 40);
	//testLineClip();
	//testPolygonClip();

	screePt p1, p2, p3;
	p1.setCoords(-100, 0);
	p2.setCoords(0, 100);
	p3.setCoords(100, -100);
	drawTrangle(p1, p2, p3);

	glFlush();
}

//输入形绘制
void inputDisplayFunc(void)
{
	GLuint curveNum;

	glClear(GL_COLOR_BUFFER_BIT);
	cout << "\n Enter the integer value corresponding to \n";
	cout << "one of the following curve names.\n";
	cout << "Press any other key to exit.\n";
	cout << "\n 1.limacon, 2-cardioid, 3-threeLeaf, 4-fourLeaf, 5-spiral:";
	cin >> curveNum;

	if (glIsList(curveNum) == GL_TRUE)
		glCallList(curveNum);
	else
		exit(0);

	glFlush();
}
//2d转换绘制
void matrixDisplayFunc()
{
	wcPt2D p1, p2, p3;
	p1.setCoords(-50.0, -50.0);
	p2.setCoords(-150.0, -50.0);
	p3.setCoords(-100.0, 50.0);
	wcPt2D verts[3] = { p1, p2, p3 };
	GLint nVerts = 3;

	wcPt2D centriodPt;

	GLint k, xSum = 0, ySum = 0;

	for (k = 0; k < nVerts; ++k)
	{
		xSum += verts[k].x;
		ySum += verts[k].y;
	}

	centriodPt.x = GLfloat(xSum) / GLfloat(nVerts);
	centriodPt.y = GLfloat(ySum) / GLfloat(nVerts);

	wcPt2D pivPt, fixedPt;
	pivPt = centriodPt;
	fixedPt = centriodPt;

	GLfloat tx = 200.0, ty = 0.0;
	GLfloat sx = 1, sy = 1;
	GLdouble theta = PI / 2.0;

	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.0f, 1.0f);//设定显示对象的颜色
	triangle2D(verts, nVerts);

	matrix3x3SetIdentity(matComposite);
	scale2D(sx, sy, fixedPt);
	rotate2D(pivPt, theta);
	translate2D(tx, ty);

	transformVerts2D(nVerts, verts);

	glColor3f(1.0, .0, .0);
	triangle2D(verts, nVerts);

	glFlush();
}
//3d转换绘制
float angle = 0;
void matrixDisplay3DFunc()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.1f, 0.1f, 0.1f);

	lineBres(-400, 0, 400, 0);
	lineBres(0, -300, 0, 300);

	int nVerts = 3;
	wcPt3D p1, p2, p3;
	p1.setCoords(-50, 0.0, 0);
	p2.setCoords(0.0, 100.0, 0);
	p3.setCoords(50.0, 50.0, 0);
	wcPt3D verts[3] = { p1, p2, p3 };

	//glColor3f(0.0f, 0.0f, 1.0f);
	//triangle3D(verts, nVerts);

	//translate3D(100, 0, 0);

	wcPt3D viewOrigin;
	viewOrigin.setCoords(0, 0, 0);
	wcPt3D viewLookAt;
	viewLookAt.setCoords(0, 50, 0);

	matrix4x4SetIdentity(matComposite3D);
	rotate3D(viewOrigin, viewLookAt, PI);
	transformVerts3D(nVerts, verts);
	
	//identity Camera Matrix
	wcPt3D cameraPos, cameraLookAt, cameraUpDir;
	cameraPos.setCoords(200*sin(angle), 0, 200*cos(angle));
	cameraLookAt.setCoords(0, 0, 0);
	cameraUpDir.setCoords(0, 1, 0);
	Matrix4x4 cameraMatrix;
	generateCameraModel(cameraMatrix, cameraPos, cameraLookAt, cameraUpDir);

	matrix4x4SetIdentity(matComposite3D);
	matrix4x4PreMultiply(cameraMatrix, matComposite3D);
	transformVerts3D(nVerts, verts);
	
	//identity Project Matrix
	float fov = PI*3/9, aspect = 1;
	float zNear = 0, zFar = 100;
	Matrix4x4 projectMatrix;
	generateProjectModel(projectMatrix, fov, aspect, zNear, zFar);

	matrix4x4SetIdentity(matComposite3D);
	matrix4x4PreMultiply(projectMatrix, matComposite3D);
	transformVerts3D(nVerts, verts);
	
	//homoneous transformation
	int k;
	for (k = 0; k < nVerts; ++k)
		verts[k].homogeneous();

	//identity Screen Matrix
	wcPt3D screenCenter;
	screenCenter.setCoords(0, 0, 0);
	Matrix4x4 screenMatrix;
	generateScreenModel(screenMatrix, screenCenter, 400, 400);

	matrix4x4SetIdentity(matComposite3D);
	matrix4x4PreMultiply(screenMatrix, matComposite3D);
	transformVerts3D(nVerts, verts);
	
	glColor3f(1.0f, 0.0f, 0.0f);
	triangle3D(verts, nVerts);

	glFlush();
}

void idleAnimFunc()
{
	angle += PI/3600;
	glutPostRedisplay();
}

void init(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);//设定背景颜色

	glMatrixMode(GL_PROJECTION);//设定投影模式
	glLoadIdentity();
	gluOrtho2D(-100.0, 100.0, -100.0, 100.0);
}

void winReshapeFunc(int newWidth, int newHeight)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(newWidth*-0.5, (GLdouble)newWidth * 0.5, (GLdouble)newHeight*-0.5, (GLdouble)newHeight*0.5);

	glClear(GL_COLOR_BUFFER_BIT);
}

int main(int argc, char** argv) 
{
	glutInit(&argc, argv);//初始化

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//设定窗口的缓存和颜色模型, GLUT_RGB-颜色编码，GLUT_INDEX-颜色表

	glutInitWindowPosition(100, 100);//设定窗口位置

	glutInitWindowSize(800, 600);//设定窗口大小

	glutCreateWindow("Computer Graphic Explore");//设置窗口标题

	init();
	initDrawList();
	//glutDisplayFunc(rasterFunc);
	//glutDisplayFunc(drawFunc);
	//glutDisplayFunc(inputDisplayFunc);//指定显示内容
	//glutDisplayFunc(matrixDisplayFunc);
	glutDisplayFunc(matrixDisplay3DFunc);
	glutIdleFunc(idleAnimFunc);
	glutReshapeFunc(winReshapeFunc);//显示窗口重定形
	glutMainLoop();
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
