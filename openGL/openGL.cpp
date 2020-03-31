// openGL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
//#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

using namespace std;

//#pragma comment(lib, "glew32.lib")

const double TWO_PI = 6.2831853;
const double PI = 3.1415926;

GLsizei winWidth = 600, winHeight = 500;
GLint xRaster = 25, yRaster = 100;
GLubyte label[36] = { 'J', 'a', 'n',		'F', 'e', 'b',		'M', 'a', 'r',
					  'A', 'p', 'r',		'M', 'a', 'y',		'J', 'u', 'n',
					  'J', 'u', 'l',		'A', 'u', 'g',		'S', 'e', 'p',
					  'O', 'c', 't',		'N', 'o', 'v',		'D', 'e', 'c',
};

GLint dataValue[12] = { 420, 342, 324, 310, 262, 185, 190, 196, 217, 240, 312, 438 };

GLfloat xwcMin = 0.0, xwcMax = 225.0;
GLfloat ywcMin = 0.0, ywcMax = 225.0;

typedef GLfloat Matrix3x3[3][3];
Matrix3x3 matComposite;

class screePt
{
private:
	GLint x;
	GLint y;

public:
	screePt()
	{
		x = y = 0;
	}

	void setCoords(GLint xValue, GLint yValue)
	{
		x = xValue;
		y = yValue;
	}

	GLint getX() const
	{
		return x;
	}

	GLint getY()const
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
	GLfloat x, y;
	wcPt2D() { x = y = 0.0; }
	void setCoords(GLfloat vx, GLfloat vy) { x = vx; y = vy; }
	~wcPt2D() 
	{ 
		cout << "destroy wcpt2D'x:" << x << " y:" << y << endl;
	}
};

inline int roundPlus(const float a) { return int(a + 0.5f); }

void setPixel(GLint x, GLint y)
{
	glBegin(GL_POINTS);
	glVertex2i(x, y);
	glEnd();
}

void triangle(wcPt2D* verts, GLint nVerts)
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

void swapPts(wcPt2D* p1, wcPt2D* p2)
{
	wcPt2D temp;
	temp = *p1;
	*p1 = *p2;
	*p2 = temp;
}

#pragma region Matrix3x3

//
void matrix3x3SetIdentity(Matrix3x3 matIdent3x3)
{
	GLint row, col;
	for (row = 0; row < 3; ++row)
	{
		for (col = 0; col < 3; ++col)
		{
			matIdent3x3[row][col] = (row == col);
		}
	}
}

void matrix3x3PreMultiply(Matrix3x3 m1, Matrix3x3 m2)
{
	GLint row, col;
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

void translate2D(GLfloat tx, GLfloat ty)
{
	Matrix3x3 matTrans;
	matrix3x3SetIdentity(matTrans);

	matTrans[0][2] = tx;
	matTrans[1][2] = ty;

	matrix3x3PreMultiply(matTrans, matComposite);
}

void rotate2D(wcPt2D pivotPt, GLfloat theta)
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

void scale2D(GLfloat sx, GLfloat sy, wcPt2D fixedPt)
{
	Matrix3x3 matScale;
	matrix3x3SetIdentity(matScale);

	matScale[0][0] = sx;
	matScale[0][2] = (1 - sx) * fixedPt.x;

	matScale[1][1] = sy;
	matScale[1][2] = (1 - sy) * fixedPt.y;

	matrix3x3PreMultiply(matScale, matComposite);
}

void transformVerts2D(GLint nVerts, wcPt2D* verts)
{
	GLint k;
	GLfloat temp;

	for (k = 0; k < nVerts; ++k)
	{
		temp = matComposite[0][0] * verts[k].x + matComposite[0][1] * verts[k].y + matComposite[0][2];
		verts[k].y = matComposite[1][0] * verts[k].x + matComposite[1][1] * verts[k].y + matComposite[1][2];
		verts[k].x = temp;
	}
}

#pragma endregion

#pragma region Draw Sampler

void lineSegment(screePt p1, screePt p2)
{
	glBegin(GL_LINES);
	glVertex2i(p1.getX(), p1.getY());
	glVertex2i(p2.getX(), p2.getY());
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

GLuint drawCurve(GLint curveNum)
{
	GLuint regHex = glGenLists(1);
	glNewList(regHex, GL_COMPILE);

	const GLint a = 175, b = 60;

	GLfloat r, theta, dtheta = 1.0 / (float)(a);
	GLint x0 = 200, y0 = 250;

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

void lineDDA(int x0, int y0, int xEnd, int yEnd)
{
	int dx = xEnd - x0;
	int dy = yEnd - y0;
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
	int dx = abs(xEnd - x0);
	int dy = abs(yEnd - y0);
	int dy2 = 2 * dy;
	int a = 2 * (dy - dx);
	int p = dy2 - dx;
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

//中点法画圆
void circlePlotPoints(GLint x, GLint y, screePt point)
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
void circleMidPoint(GLint x, GLint y, GLint radius)
{
	screePt circPt;
	GLint p = 1 - radius;
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
void ellipsePlotPoints(GLint xCenter, GLint yCenter, GLint x, GLint y)
{
	setPixel(xCenter + x, yCenter + y);
	setPixel(xCenter + x, yCenter - y);
	setPixel(xCenter - x, yCenter + y);
	setPixel(xCenter - x, yCenter - y);
}
void ellipseMidPoint(GLint xCenter, GLint yCenter, GLint xRadius, GLint yRadius)
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
void midPointGravity(GLint vx0, GLint vy0, GLint x0, GLint y0)
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

#pragma endregion

#pragma region 2D Clip

const GLint winLeftBitCode = 0x1;
const GLint winRightBitCode = 0x2;
const GLint winBottomBitCode = 0x4;
const GLint winTopBitCode = 0x8;

inline GLint inside(GLint code) { return GLint(!code); }//取反
inline GLint reject(GLint code1, GLint code2) { return GLint(code1 & code2); }//且
inline GLint accept(GLint code1, GLint code2) { return GLint(!(code1 | code2)); }//异或

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
	GLint done = false, plotLine = false;

	GLfloat m;

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
							if (p2.x != p1.x)
								p1.x += (winMin.y - p1.y) / m;
							p1.y = winMin.y;
						}
						else
						{
							if (code1 & winTopBitCode)
							{
								if (p2.x != p1.x)
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

GLint clipTest(GLfloat p, GLfloat q, GLfloat* uIn, GLfloat* uOut)
{
	GLfloat r;
	GLint returnValue = true;

	if (p < .0)//入边
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
		if (p > .0)//出边
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
		else//与边界平等
			if (q < .0)
				returnValue = false;
	}
	return returnValue;
}

void lineClipLiangBarsk(wcPt2D& winMin, wcPt2D& winMax, wcPt2D& p1, wcPt2D& p2)
{
	GLfloat uIn = 0, uOut = 1, dx = p2.x - p1.x, dy;
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
	wcPt2D winMin;
	winMin.x = 0; winMin.y = 0;
	wcPt2D winMax;
	winMax.x = 100; winMax.y = 100;

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
const GLint nClip = 4;

//是否在边界区域内
GLint inside(wcPt2D& p, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax)
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
GLint cross(wcPt2D& p1, wcPt2D& p2, Boundary edeg, wcPt2D& winMin, wcPt2D& winMax)
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
	GLfloat m;

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
		first[edeg] = &p;
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

GLint polygonClipSuthHodg(wcPt2D& winMin, wcPt2D& winMax, GLint n, wcPt2D* pIn, wcPt2D* pOut)
{
	wcPt2D* first[nClip] = { 0, 0, 0, 0 };
	wcPt2D last[nClip];
	GLint k, cnt = 0;

	for (k = 0; k < n; ++k)
		clipPoint(pIn[k], Left, winMin, winMax, pOut, &cnt, first, last);

	wcPt2D pt;
	for (k = 0; k < nClip; ++k)
	{
		pt.setCoords(first[k]->x, first[k]->y);
		cout << "first [" << k << "]" << " Point x: " << pt.x << endl;
	}
	closeClip(winMin, winMax, pOut, &cnt, first, last);
	return cnt;
}

void testPolygonClip()
{
	wcPt2D winMin, winMax;
	winMin.x = 0; winMin.y = 0;
	winMax.x = 100; winMax.y = 100;
	lineDDA(0, 0, 100, 0);
	lineDDA(100, 0, 100, 100);
	lineDDA(0, 100, 100, 100);
	lineDDA(0, 0, 0, 100);

	GLint nVerts = 3;
	wcPt2D p1, p2, p3;
	p1.setCoords(-25.0, 10.0);
	p2.setCoords(125.0, 10.0);
	p3.setCoords(50.0, 125.0);
	wcPt2D verts[3] = { p1, p2, p3};

	wcPt2D p4, p5, p6, p7, p8, p9, p10;
	wcPt2D outVerts[7] = { p4, p5, p6, p7, p8, p9, p10 };

	GLint outCnt = polygonClipSuthHodg(winMin, winMax, nVerts, verts, outVerts) - 1;
	triangle(outVerts, outCnt);
}

#pragma endregion

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
	//drawColorTriangles();
	//lineBres(20, 10, 330, 218);
	//circleMidPoint(200, 200, 100);
	//ellipseMidPoint(200, 200, 100, 50);
	//midPointGravity(200, 800, 40, 40);
	//testLineClip();
	testPolygonClip();

	glFlush();
}


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
	triangle(verts, nVerts);

	matrix3x3SetIdentity(matComposite);
	scale2D(sx, sy, fixedPt);
	rotate2D(pivPt, theta);
	translate2D(tx, ty);

	transformVerts2D(nVerts, verts);

	glColor3f(1.0, .0, .0);
	triangle(verts, nVerts);

	glFlush();
}

void init(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);//设定背景颜色

	glMatrixMode(GL_PROJECTION);//设定投影模式
	glLoadIdentity();
	gluOrtho2D(-200.0, 200.0, -200.0, 200.0);
}

void winReshapeFunc(int newWidth, int newHeight)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(newWidth*-0.5, (GLdouble)newWidth * 0.5, (GLdouble)newHeight*-0.5, (GLdouble)newHeight*0.5);

	glClear(GL_COLOR_BUFFER_BIT);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);//初始化

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//设定窗口的缓存和颜色模型, GLUT_RGB-颜色编码，GLUT_INDEX-颜色表

	glutInitWindowPosition(100, 100);//设定窗口位置

	glutInitWindowSize(400, 400);//设定窗口大小

	glutCreateWindow("OpenGL Explore");//设置窗口标题

	init();
	initDrawList();
	glutDisplayFunc(drawFunc);
	//glutDisplayFunc(inputDisplayFunc);//指定显示内容
	//glutDisplayFunc(matrixDisplayFunc);//
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
