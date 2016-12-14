/*
* hough.h
*
*  Created on: 15. nov. 2016
*      Author: Gasper
*/

#ifndef DEBUG_HOUGH_H_
#define DEBUG_HOUGH_H_



#endif /* DEBUG_HOUGH_H_ */

double PI = 3.14159265;
double outOfBounds = 1000000;
double box[3][3] = { { 1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0 }, { 1.0 / 8.0, 1.0 / 4.0, 1.0 / 8.0 }, { 1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0 } };


struct FuncArgs{
	int threadStartRow;
	int threadEndRow;
	int tfRows;
	int tfCols;
	int orgRows;
	int orgCols;
	int numLines;
	double diag;
	int numThetas;
	int minNumPixelsOnLine;
	double** arrIn;
	double**arrOut;
	double** lines;

};


//arithmetic-logic functions
void* houghField(void* arg);
void fixEdges(FuncArgs *args);
void blur(void *args);
void keepMaximas(FuncArgs * args);
void getMaximas(FuncArgs * args);
void getDesiredLines(FuncArgs* args);
double** getLines(FuncArgs * args);
void sortLines(FuncArgs * args);

void mergeLines(FuncArgs* args);
void mergeTwoLines(double* a, double* b);
void rotateLines(FuncArgs *args);
double** removeUndesiredLines(FuncArgs *args);
double** removeOutOfFieldLines(FuncArgs *args);
double** calculateNormalizedDistances(double** vertices);
double** removeOutOfBorderIntersections(FuncArgs* args, double** isects);
double* findIntersect(double* line0, double* line1);
double** generateIntersectTable(double ** lines, int numLines);
double calculateAngle(double theta);
double calculateX0Intersection(double theta, double r);
double* calculateLine(double theta, double r);
/*
//I/O functions
void outArr(int** field, int x, int y);
void outImgI(int** field, int x, int y);
*/
void outLines(FuncArgs *args);
void outImgD(double** field, int x, int y);
void outIntersectTable(double** table);
void outDistances(double** distances);

//shifters&helpers
void edgify(FuncArgs *args);
void polar2cartesian(double**lines);
double** makeDoublearray(int x, int y);
int** makeIntarray(int x, int y);
void zeroArray(double **arr, int x, int y);
void get2dFrom1d(int**  img, uchar* a, int rows, int cols);

double absolute(double x);

void edgify(FuncArgs *args){
	double ** img = args->arrIn;
	int i = 1;
	int rows = args->tfRows;
	int cols = args->tfCols;
	for (; i < rows - 1; i++) {
		int j = 1;
		for (; j < cols; j++) {
			if (img[i - 1][j - 1]>0 && img[i - 1][j]>0 && img[i - 1][j + 1]>0 && img[i][j - 1]>0 && img[i][j]>0 && img[i][j + 1]>0 && img[i + 1][j - 1]>0 && img[i + 1][j]>0 && img[i + 1][j]>0)
				img[i][j] = 300;
		}
	}
	i = 1;
	for (; i < rows - 1; i++) {
		int j = 1;
		for (; j < cols; j++) {
			if (img[i][j] == 300)img[i][j] = 0;
		}
	}
	args->arrOut = img;
}

void* houghField(void *arg) {
	FuncArgs* args = ((FuncArgs*)arg);
	double ** img = args->arrIn;
	double ** htr = args->arrOut;
	int rows = args->orgRows;
	int cols = args->orgCols;
	int start = args->threadStartRow;;
	int end = args->threadEndRow;

	int i = start;
	for (; i < end; i++) {
		int j = 0;
		for (; j < cols; j++) {
			if (img[i][j] < 250)
				continue;
			int theta = 0;
			for (; theta < args->numThetas; theta++) {
				double angle = PI * (((double)theta) / ((double)args->numThetas));
				double r = ((double)(i + 1)) * cos(angle) + ((double)(j + 1)) * sin(angle);
				r += args->diag;
				double diff = r - (int)r;
				htr[1 + (int)r][1 + theta] += 1 - diff;
				htr[2 + (int)r][1 + theta] += diff;
			}
		}
	}
	return arg;
}

void fixEdges(FuncArgs *args) {
	double ** img = args->arrIn;
	int rows = args->tfRows;
	int cols = args->tfCols;
	int start = args->threadStartRow;
	int end = args->threadEndRow;

	int i = 0;
	for (; i < rows; i++) {
		img[i][0] = img[i][1];
		img[i][cols - 1] = img[i][cols - 2];
	}
	i = 0;
	for (; i < cols; i++) {
		img[0][i] = img[1][i];
		img[rows - 1][i] = img[rows - 2][i];
	}
	args->arrOut = args->arrIn;
}

double** calculateNormalizedDistances(double** v) {
	int n = (int)v[0][0];
	int numOfEdges = n * (n - 1) / 2;
	double* distances = new double[numOfEdges];
	int i = 1, edgeIter = 0;
	for (; i <= n; i++) {
		int j = 1;
		for (; j < i; j++) {
			distances[edgeIter++] = sqrt(
				(v[i][0] - v[j][0]) * (v[i][0] - v[j][0])
				+ (v[i][1] - v[j][1]) * (v[i][1] - v[j][1]));
		}
	}

	double**normalizedDistances = new double*[numOfEdges + 1];
	normalizedDistances[0] = new double[1];
	normalizedDistances[0][0] = numOfEdges + 1;
	i = 1;
	for (; i <= numOfEdges; i++)
		normalizedDistances[i] = new double[numOfEdges];
	i = 1;
	for (; i <= numOfEdges; i++) {
		double maxVal = distances[i - 1];
		int j = 0;
		for (; j < numOfEdges; j++) {
			normalizedDistances[i][j] = distances[j] / maxVal;
		}
	}
	return normalizedDistances;
}

double** generateIntersectTable(FuncArgs* args) {
	double** lines = args->lines;
	int numLines = args->numLines;
	double** iTable = new double*[1];
	iTable[0] =new double[1];
	iTable[0][0] = 0;
	int i = 0;
	for (; i < numLines; i++) {
		int j = 0;
		for (; j < i; j++) {
			int iTabIter = iTable[0][0] + 1;

			double** newItab = new double*[iTabIter + 1];
			for (int m=0; m < iTabIter; m++)
				newItab[m] = iTable[m];
			iTable = newItab;

			iTable[iTabIter] = findIntersect(lines[i], lines[j]);
			iTable[0][0]++;
		}
	}
	return removeOutOfBorderIntersections(args, iTable);
}

double* findIntersect(double* line0, double* line1) {
	double* vals = new double[2];

	if (line0[0] != outOfBounds && line1[0] != outOfBounds) {
		double xintersect = (line1[1] - line0[1]) / (line0[0] - line1[0]);
		double yintersect = line1[1] + xintersect * line1[0];
		vals[0] = xintersect;
		vals[1] = yintersect;
	}
	else if (line1[0] != outOfBounds) {
		vals[0] = line0[1];
		vals[1] = line1[1] + line0[1] * line1[0];
	}
	else if (line0[0] != outOfBounds) {
		vals[0] = line1[1];
		vals[1] = line0[1] + line1[1] * line0[0];
	}
	else {
		vals[0] = outOfBounds;
		vals[1] = outOfBounds;
	}
	return vals;
}

double** removeOutOfBorderIntersections(FuncArgs*args, double** isects){
	int i = 1, allISects = isects[0][0];
	for (; i <= isects[0][0]; i++) {
		if (isects[i][0]<-10 || isects[i][0]> args->orgRows || isects[i][1]<0 || isects[i][1]>args->orgCols){
			allISects--;
			isects[i][0] = -1;
			isects[i][1] = -1;
		}
	}
	double** visibleISects = new double*[allISects + 1];
	visibleISects[0] = new double[1];
	visibleISects[0][0] = allISects;
	int visIter = 1;
	i = 1;
	for (; i <= isects[0][0]; i++){
		if (isects[i][0] == -1 && isects[i][1] == -1)
			continue;
		visibleISects[visIter++] = isects[i];
	}
	return visibleISects;
}

void getDesiredLines(FuncArgs* args){

	getLines(args);
	sortLines(args);
	args->lines = removeUndesiredLines(args);
}

double** getLines(FuncArgs * args) {
	//potential lines are visible in the houghTransform field because they are the maximal value there
	double ** htr = args->arrIn;
	int numThetas = args->numThetas;
	int diag = args->diag;
	//PRAV
	getMaximas(args);
	double** maximas = args->lines;
	int numMaximas = args->numLines;
	args->lines = maximas;
	keepMaximas(args);
	getMaximas(args); 
	maximas = args->lines;
	numMaximas = args->numLines;
	int i = 0;
	for (; i < numMaximas; i++) {
		double theta = maximas[i][0] - 1;
		double r = maximas[i][1] - diag - 1;
		maximas[i][0] = theta / (double)numThetas;
		maximas[i][1] = r;
	}
	args->lines = maximas;
	return maximas;
}

void getMaximas(FuncArgs * args){
	double ** htr = args->arrIn;
	double ** ret = args->arrOut;
	int rows = args->tfRows;
	int cols = args->tfCols;
	int start = args->threadStartRow;
	int end = args->threadEndRow;
	int lineCount = 0;
	double** lines = new double*[0];
	int i = start;
	for (; i < end; i++) {
		int j = 0;
		for (; j < cols; j++) {
			if (htr[i][j] <= (double)args->minNumPixelsOnLine)
				continue;
			if (htr[i][j] > htr[i - 1][j - 1] && htr[i][j] > htr[i - 1][j]
				&& htr[i][j] > htr[i - 1][j + 1]
				&& htr[i][j] > htr[i][j - 1] && htr[i][j] >= htr[i][j + 1]
				&& htr[i][j] >= htr[i + 1][j - 1]
				&& htr[i][j] >= htr[i + 1][j]
				&& htr[i][j] >= htr[i + 1][j + 1]) {

				double** newLines = new double*[lineCount + 1];
				for (int m = 0; m < lineCount; m++)
					newLines[m] = lines[m];

				free(lines);
				lines = newLines;
				lines[lineCount] = new double[3];
				lines[lineCount][0] = j;
				lines[lineCount][1] = i;
				lines[lineCount][2] = htr[i][j];
				lineCount++;
			}
		}
	}
	args->lines = lines;
	args->numLines = lineCount;
}

void keepMaximas(FuncArgs * args){
	double** blur = args->arrIn;
	double** htr = args->arrOut;
	double** maximas = args->lines;
	int numMaximas = args->numLines;
	int rows = args->tfRows;
	int cols = args->tfCols;

	int l = 0;
	for (; l<numMaximas; l++){
		int j = (int)maximas[l][0];
		int i = (int)maximas[l][1];
		int k = -1;
		for (; k <= 1; k++){
			int newi = i + k;
			if (newi<0 || newi >= rows)
				continue;
			int l = -1;
			for (; l <= 1; l++){
				int newj = j + l;
				if (newj<0 || newj >= cols)
					continue;
				blur[newi][newj] = -1;
			}
		}
	}
	int i = 0;
	for (; i<rows; i++){
		int j = 0;
		for (; j<cols; j++){
			if (blur[i][j] == -1)
				blur[i][j] = htr[i][j];
			else
				blur[i][j] = 0;
		}
	}
}

void outImgD(double** field, int x, int y) {
	double max = 0;
	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			double f = field[i][j];
			if (f > max)
				max = field[i][j];
		}
	}
	i = 0;
	printf("%d %d \n", x, y);
	printf("\n");
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			printf("%d ", (int)(255.0 * field[i][j] / max));
			//printf("%f ",field[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
}


void polar2cartesian(FuncArgs *args) {
	double** lines = args->lines;
	int numLines = args->numLines;
	int i = 0;
	for (; i < numLines; i++) {
		if (absolute(lines[i][0] - 0.5) < 0.01 || absolute(lines[i][0] + 0.5) < 0.01) {
			//case of vertical line(cartesian angle for vertical line is infinity)
			lines[i][0] = outOfBounds;
		}
		else {
			double* cartesianLine = calculateLine(lines[i][0], lines[i][1]);
			free(lines[i]);
			lines[i] = cartesianLine;
		}

	}
}

double* calculateLine(double theta, double r) { //y= ax+b => ret[0]=a;ret[1]=b
	double* vals = new double[2];
	vals[0] = calculateAngle(theta);
	vals[1] = calculateX0Intersection(theta, r);
	return vals;
}

double calculateAngle(double theta) {
	double angle = -tan(PI * theta);
	return angle;
}

double calculateX0Intersection(double theta, double r) {
	double cosAng = cos(PI * theta);
	double intersect = r / cosAng;
	return intersect;
}


void blur(void * arg) {
	FuncArgs* args = ((FuncArgs*)arg);

	double ** img = args->arrIn;
	double ** ret = args->arrOut;
	int rows = args->tfRows;
	int cols = args->tfCols;
	int start = args->threadStartRow;
	int end = args->threadEndRow;
	int i = start;
	for (; i < end - 1; i++) {
		if (i == 0 || i == rows - 1)continue;
		int j = 1;
		for (; j < cols - 1; j++) {
			if (j == 1 || j == cols - 1)continue;
			ret[i][j] = 0;
			ret[i][j] = img[i - 1][j - 1] * box[0][0];
			ret[i][j] += img[i - 1][j] * box[0][1];
			ret[i][j] += img[i - 1][j + 1] * box[0][2];
			ret[i][j] += img[i][j - 1] * box[1][0];
			ret[i][j] += img[i][j] * box[1][1];
			ret[i][j] += img[i][j + 1] * box[1][2];
			ret[i][j] += img[i + 1][j - 1] * box[2][0];
			ret[i][j] += img[i + 1][j] * box[2][1];
			ret[i][j] += img[i + 1][j + 1] * box[2][2];
		}
	}
}


void sortLines(FuncArgs * args) {
	double** vals = args->lines;
	int numOfLines = args->numLines;

	int i = 0;
	for (; i < numOfLines; i++) {
		int j = i + 1;
		for (; j < numOfLines; j++) {
			if (vals[i][2] < vals[j][2]) {
				double *tmp = vals[i];
				vals[i] = vals[j];
				vals[j] = tmp;
			}
		}
	}
}

void rotateLines(FuncArgs * args) {
	double** lines = args->lines;
	double numOfLines = args->numLines;
	int i = 0;
	for (; i < numOfLines; i++) {
		if (lines[i][0] > 0.5 && lines[i][1] < 0) {
			lines[i][0] -= 1.0;
			lines[i][1] *= -1.0;
		}
	}
}

void mergeLines(FuncArgs* args) {
	double ** lines = args->lines;
	double numOfLines = args->numLines;
	double acceptableAngleDiff = 0.2;
	double acceptableRDiff = 0.08*args->diag;
	int i = 0;
	for (; i < numOfLines; i++) {
		if (lines[i][0] == 0 && lines[i][1] == -1) {
			continue;
		}
		int j = i + 1;
		for (; j < numOfLines; j++) {
			if (lines[j][0] == 0 && lines[j][1] == -1)
				continue;
			double absA = absolute(lines[i][0] - lines[j][0]);
			double absB = absolute(lines[i][1] - lines[j][1]);

			if (absA < acceptableAngleDiff && absB <= acceptableRDiff) {
				mergeTwoLines(lines[i], lines[j]);	//are simmilar!!! Merge them in to one line!
				i = -1;
				break;
			}
		}
	}
	double** cleanedLines = removeOutOfFieldLines(args);
	args->lines = cleanedLines;
}

void mergeTwoLines(double* a, double* b) {
	//printf("MERGING LINES:\n1:%f;%f\n2:%f;%f\n",a[0],a[1],b[0],b[1]);
	double sumIntensity = a[2] + b[2];
	double aFactor = a[2] / sumIntensity;
	double bFactor = b[2] / sumIntensity;
	a[0] = a[0] * aFactor + b[0] * bFactor;
	a[1] = a[1] * aFactor + b[1] * bFactor;
	a[2] = sqrt(a[2] * a[2] + b[2] * b[2]);
	b[0] = 0;
	b[1] = -1;
	b[2] = 0;
	//printf("m:%f;%f\n\n",a[0],a[1]);
	//fflush(stdout);
}

double** removeOutOfFieldLines(FuncArgs *args) {
	double** lines = args->lines;
	int numOfLines = args->numLines;

	int i = 0, newNumOfLines = 0;
	for (; i < numOfLines; i++) {
		if (lines[i][0] != 0 && lines[i][0] != -1 && lines[i][2] != 0)
			newNumOfLines++;
	}
	double** linesNew = new double*[newNumOfLines];
	i = 0;
	int iNew = 0;
	for (; i < numOfLines; i++) {
		if (lines[i][0] != 0 && lines[i][0] != -1 && lines[i][2] != 0)
			linesNew[iNew++] = lines[i];
		else
			free(lines[i]);
	}
	free(lines);
	args->numLines = newNumOfLines;
	args->lines = linesNew;
	return linesNew;
}

double** removeUndesiredLines(FuncArgs *args) {
	double** lines = args->lines;

	int numOfLines = args->numLines;
	int sum = 0.0, max = 0.0;
	int i = 1;
	for (; i < numOfLines; i++) {
		sum += lines[i][2];
		if (lines[i][2]>max)max = lines[i][2];
	}
	//double treshold = sum / (numOfLines - 1.0);
	double treshold = max*0.4;
	i = 1;
	for (; i < numOfLines; i++) {
		if (lines[i][2] < treshold) {
			lines[i][0] = 0;
			lines[i][1] = -1;
			lines[i][2] = 0;
		}
	}
	return removeOutOfFieldLines(args);
}

double** makeDoublearray(int x, int y) {
	double **ret = new double*[x];
	int i = 0;
	for (; i < x; i++) {
		ret[i] = new double[y];
		int j = 0;
		for (; j < y; j++) {
			ret[i][j] = 0;
		}
	}
	return ret;
}

int** makeIntarray(int x, int y) {
	int **ret = new int*[x];
	int i = 0;
	for (; i < x; i++) {
		ret[i] = new int[y];
		int j = 0;
		for (; j < y; j++) {
			ret[i][j] = 0;
		}
	}
	return ret;
}

void zeroArray(double **arr, int x, int y) {
	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			arr[i][j] = 0;
		}
	}
}

void get2dFrom1d(int** table, uchar* a, int rows, int cols){
	int i=0;
	for(;i<rows;i++){
		int j=0;
		for(;j<cols;j++){
			table[i][j]=(int)a[(j*3)+(i*cols*3)];
		}
	}
}

void outLines(FuncArgs *args) {
	//printf("FORMAT OF OUTPUT:\n [i]: A:[a], B:[b], V:[v]; \n	where i is the number of the line coordinates,\n	Form: y=ax+b \n	and v denotes the number of pixels on that line\n\n");
	int numLines = args->numLines;
	double** lines = args->lines;
	int i = 0;
	printf("Number of lines:%d\n", (int)numLines);
	for (; i < numLines; i++) {
		printf("%d: A:%f, B:%f, V:%f\n", i, lines[i][0], lines[i][1],
			lines[i][2]);
	}
	printf("\n");
	fflush(stdout);
}



void outIntersectTable(double** table) {
	printf("Calculated vertices of figure:\n");
	double numIntersects = table[0][0] + 1;
	int i = 1;
	for (; i < numIntersects; i++) {
		printf("[%d, %d], \n", (int)table[i][0], (int)table[i][1]);
	}
	printf("\n");

}

void outDistances(double** distances) {
	int numOfDists = distances[0][0];
	int i = 1;
	printf("Calculated distances of figure:\n");
	for (; i < numOfDists; i++) {
		int j = 0;
		for (; j < numOfDists - 1; j++) {
			printf("%f,	", distances[i][j]);
		}
		printf("\n");
	}
}

double absolute(double x){
	if(x<0)return -x;
	return x;
}


/*
void outArr(int** field, int x, int y) {
int i = 0;
for (; i < x; i++) {
int j = 0;
for (; j < y; j++) {
printf("%d ", field[i][j]);
}
printf("\n");
}
}

void outImgI(int** field, int x, int y) {
double max = 0;
int i = 0;
for (; i < x; i++) {
int j = 0;
for (; j < y; j++) {
if (field[i][j] > max)
max = field[i][j];
}
}

printf("int[][] img={");
i = 0;
for (; i < x; i++) {
int j = 0;
printf("{");
for (; j < y; j++) {
printf("%d, ", (int) (255 * ((double) field[i][j])));
}
printf("},\n");
}
printf("};\n");
}

*/
