#include <cmath>
#include <iostream>
using std::FILE;

#define PI  3.14159265359
//==================解算EO// Exterior Orientation 函数==============================//
//      函数名称:  ExtractEO
//		功能：根据POS文件解算EO
//      参数：pPOSFile：	POS文件全路径
//			  nLines：		POS行数
//			  THETA:		坐标 
//      返回：long
//      说明：
//===============================================================//
long  ExtractEO(char *pPOSFile, int nLines, char *pEOFile, THREEDPOINT THETA, float *pfRoll, float *pfHeight, float *pfHead, float *fSetupAngle)
{
	long   IError = 0;
	double  *pfPosData = new double[6 * nLines + 100];
	double  dB = 0, dL = 0, dH = 0;
	int nQuarFlag = 0;
	FILE *fEO = NULL;



	THREEDPOINT ANGLETHETA, POSISTTHETA;

	ANGLETHETA.dX = fSetupAngle[0]; ANGLETHETA.dY = fSetupAngle[1]; ANGLETHETA.dZ = fSetupAngle[2];
	POSISTTHETA.dX = THETA.dX;      POSISTTHETA.dY = THETA.dY;      POSISTTHETA.dZ = THETA.dZ;

	IError = ReadPartPOS(pPOSFile, nLines, dB, dL, dH, pfPosData);//从POS文件中读取所需的参数，计算经纬度高度的均值
	if (IError>0)
	{
		goto ErrEnd;
	}

	double EMMatrix[9];
	EMMatrix[0] = -sin(dL);
	EMMatrix[1] = cos(dL);
	EMMatrix[2] = 0;

	EMMatrix[3] = -sin(dB)*cos(dL);
	EMMatrix[4] = -sin(dB)*sin(dL);
	EMMatrix[5] = cos(dB);

	EMMatrix[6] = cos(dB)*cos(dL);
	EMMatrix[7] = cos(dB)*sin(dL);
	EMMatrix[8] = sin(dB);

	//物方坐标系以测区的重心为原点, 方向为东(X)北(Y)天(Z)的局部切平面(大地高为0)坐标系
	THREEDPOINT XYZPnt;
	BLHToXYZ(dB, dL, 0, XYZPnt);

	nQuarFlag = EOQuadrant(pfPosData, nLines, EMMatrix, XYZPnt);  //计算象限失败

	if (nQuarFlag == 0)
	{
		IError = 40310; //生成EO文件失败
		goto ErrEnd;
	}
	fEO = fopen(pEOFile, "w");
	if (fEO == NULL)
	{
		IError = 40311; //生成EO文件失败
		goto ErrEnd;
	}


	for (int i = 0; i<nLines; i++)
	{
		EOMatrixTurn(pfPosData, i, XYZPnt, nQuarFlag, EMMatrix, ANGLETHETA, POSISTTHETA, fEO);
	}


ErrEnd:
	if (fEO)
	{
		fclose(fEO);
		fEO = NULL;
	}
	delete[]pfPosData;

	return IError;
}
//========================读取POS函数==============================//
//      函数名称:  ReadPartPOS
//		功能：从POS文件中读取所需的参数，计算经纬度高度的均值
//      参数：pPOSFile：	POS文件全路径
//			  nLines：		POS行数
//			  dB、dL、dH:	经纬度、高度
//            pfPosData:    读取到的POS数据
//      返回：long
//      说明：
//===============================================================//
long  ReadPartPOS(char *pPOSFile, long nLines, double &dB, double &dL, double &dH, double *pfPosData)
{
	long   IError = 0;
	int    i;
	FILE   *fpin = NULL;
	char   cA[111], cB[111];
	double fReadData[20];
	double fTemp[6];
	char   cTempchar[2048];

	if ((fpin = fopen(pPOSFile, "r")) == NULL)
	{
		IError = 40307;  //打开POS文件失败
		goto ErrEnd;
	}

	dB = dL = dH = 0.0;

	for (i = 0; i < nLines; i++)
	{
		if (feof(fpin))
		{
			IError = 40308; //读取POS文件失败
			goto ErrEnd;
		}
		fgets(cTempchar, 2048, fpin);


		IError = sscanf(
			cTempchar,
			"%lf%lf%lf%lf%lf%s%lf%lf%lf%s%lf%lf%lf%lf%lf%lf%lf",
			&fReadData[0], &fReadData[1], &fReadData[2], &fReadData[3],
			&fReadData[4], cA, &fTemp[0], &fTemp[1], &fTemp[2], cB,
			&fTemp[3], &fTemp[4], &fTemp[5], &fReadData[7],
			&fReadData[8], &fReadData[9], &fReadData[10]);
		if (IError == 0)
		{
			IError = 40309; //解析POS文件失败
			goto ErrEnd;
		}
		fReadData[5] = fTemp[0] + fTemp[1] / double(60) + fTemp[2] / double(3600);
		fReadData[6] = fTemp[3] + fTemp[4] / double(60) + fTemp[5] / double(3600);

		if (fReadData[6] < -180.0 || fReadData[6] > 180.0)
		{
			IError = 40316;//经度有问题，不是-180到180度之间
			goto ErrEnd;
		}
		if (fReadData[5] <= -90.0 || fReadData[5] >= 90.0)
		{
			IError = 40317;//纬度有问题，不是-90到90度之间
			goto ErrEnd;
		}

		dB += fReadData[5];
		dL += fReadData[6];
		dH += fReadData[7];

		//fReadData[10] += 270.0;

		pfPosData[i * 6 + 0] = double(fReadData[5] * PI / 180.0);
		pfPosData[i * 6 + 1] = double(fReadData[6] * PI / 180.0);

		pfPosData[i * 6 + 2] = double(fReadData[7]);                //height

		//if(abs(fSetupAngle[0])<30)
		//	fReadData[8] += fSetupAngle[0];
		//fReadData[8] += 4.8;
		pfPosData[i * 6 + 3] = double(fReadData[8] * PI / 180.0);       //roll

		//if(abs(fSetupAngle[1])<30)
		//	fReadData[9] += fSetupAngle[1];
		//fReadData[9] += -1.8;
		pfPosData[i * 6 + 4] = double(fReadData[9] * PI / 180.0);       //pitch

		//if(abs(fSetupAngle[2])<30)
		//	fReadData[10] += fSetupAngle[2];
		//fReadData[10] += 0.0;
		pfPosData[i * 6 + 5] = double((fReadData[10])*PI / 180.0);    //heading
	}

	dB /= nLines; dB = dB*PI / 180.0;
	dL /= nLines; dL = dL*PI / 180.0;
	dH /= nLines;;

	IError = 0;

ErrEnd:
	if (fpin)
	{
		fclose(fpin);
		fpin = NULL;
	}
	return IError;
}

//========================计算象限函数==============================//
//      函数名称:  EOQuadrant
//		功能：根据POS数据，计算扫描行所处象限
//      参数：EMMatrix：	系数矩阵
//			  nLines：		POS行数
//			  XYZPoint:	经纬度、高度
//            pfPosData:    读取到的POS数据
//      返回：long
//      说明：
//===============================================================//
long  EOQuadrant(double *pfPosData, int nLines, double EMMatrix[], THREEDPOINT &XYZPoint)
{

	long nFlag = 0;

	THREEDPOINT Pntstart, PntEnd;
	BLHToXYZ(pfPosData[0], pfPosData[1], pfPosData[2], Pntstart);
	double dXstart = (Pntstart.dX - XYZPoint.dX)*EMMatrix[0] + (Pntstart.dY - XYZPoint.dY)*EMMatrix[1] + (Pntstart.dZ - XYZPoint.dZ)*EMMatrix[2];
	double dYstart = (Pntstart.dX - XYZPoint.dX)*EMMatrix[3] + (Pntstart.dY - XYZPoint.dY)*EMMatrix[4] + (Pntstart.dZ - XYZPoint.dZ)*EMMatrix[5];

	BLHToXYZ(pfPosData[6 * (nLines - 1)], pfPosData[6 * (nLines - 1) + 1], pfPosData[6 * (nLines - 1) + 2], PntEnd);			//WGS84
	double dXend = (PntEnd.dX - XYZPoint.dX)*EMMatrix[0] + (PntEnd.dY - XYZPoint.dY)*EMMatrix[1] + (PntEnd.dZ - XYZPoint.dZ)*EMMatrix[2];
	double dYend = (PntEnd.dX - XYZPoint.dX)*EMMatrix[3] + (PntEnd.dY - XYZPoint.dY)*EMMatrix[4] + (PntEnd.dZ - XYZPoint.dZ)*EMMatrix[5];

	double dx, dy;
	dx = dXend - dXstart;
	dy = dYend - dYstart;

	if (dx>0 && dy>0)
		nFlag = 1;
	else if (dx<0 && dy>0)
		nFlag = 2;
	else if (dx<0 && dy<0)
		nFlag = 3;
	else if (dx>0 && dy<0)
		nFlag = 4;

	return nFlag;
}
//========================计算旋转矩阵函数=========================//
//      函数名称:  EOMatrixTurn
//		功能：解算EO每行的旋转矩阵的参数
//      参数：pfPosData:    读取到的POS数据
//
//      返回：long
//      说明：
//===============================================================//
long  EOMatrixTurn(double *pfPosData, int nCurLine, THREEDPOINT &XYZPoint, int nQuadNum,
	double EMMatrix[], THREEDPOINT ANGLETHETA, THREEDPOINT POSISTTHETA, FILE *fEO)
{
	double ICMatrix[9];  //传感器坐标系到像空间坐标系旋转矩阵
	ICMatrix[0] = 0; 	ICMatrix[1] = -1;	ICMatrix[2] = 0;
	ICMatrix[3] = -1;	ICMatrix[4] = 0;	ICMatrix[5] = 0;
	ICMatrix[6] = 0;	ICMatrix[7] = 0;	ICMatrix[8] = -1;



	double dB, dL, dH;
	dB = pfPosData[6 * nCurLine + 0];
	dL = pfPosData[6 * nCurLine + 1];
	dH = pfPosData[6 * nCurLine + 2];

	double GEMatrix[9];//地理坐标系(Geo)到地心坐标系(Earth)的旋转矩阵
	GEMatrix[0] = -sin(dB)*cos(dL);
	GEMatrix[1] = -sin(dL);
	GEMatrix[2] = -cos(dB)*cos(dL);
	GEMatrix[3] = -sin(dB)*sin(dL);
	GEMatrix[4] = cos(dL);
	GEMatrix[5] = -cos(dB)*sin(dL);
	GEMatrix[6] = cos(dB);
	GEMatrix[7] = 0;
	GEMatrix[8] = -sin(dB);

	double dRoll, dPitch, dYaw;
	double IGMatrix[9];//IMU坐标系(IMU)到地理坐标系(Geo)的旋转矩阵
	dRoll = pfPosData[6 * nCurLine + 3];
	dPitch = pfPosData[6 * nCurLine + 4];
	dYaw = pfPosData[6 * nCurLine + 5];

	IGMatrix[0] = cos(dPitch)*cos(dYaw);
	IGMatrix[1] = sin(dRoll)*sin(dPitch)*cos(dYaw) - cos(dRoll)*sin(dYaw);
	IGMatrix[2] = cos(dRoll)*sin(dPitch)*cos(dYaw) + sin(dRoll)*sin(dYaw);
	IGMatrix[3] = cos(dPitch)*sin(dYaw);
	IGMatrix[4] = sin(dRoll)*sin(dPitch)*sin(dYaw) + cos(dRoll)*cos(dYaw);
	IGMatrix[5] = cos(dRoll)*sin(dPitch)*sin(dYaw) - sin(dRoll)*cos(dYaw);
	IGMatrix[6] = -sin(dPitch);
	IGMatrix[7] = sin(dRoll)*cos(dPitch);
	IGMatrix[8] = cos(dRoll)*cos(dPitch);


	double IMMatrix[9];
	double M1[9], M2[9];

	MatrixMuti(EMMatrix, GEMatrix, M1, 3, 3, 3);
	MatrixMuti(M1, IGMatrix, M2, 3, 3, 3);
	MatrixMuti(M2, ICMatrix, IMMatrix, 3, 3, 3);

	double dPhi, dOmega, dKappa;

	dPhi = asin(-IMMatrix[2]);
	dOmega = atan(-IMMatrix[5] / IMMatrix[8]);
	if (nQuadNum == 1 || nQuadNum == 2)
	{
		dKappa = atan(-IMMatrix[1] / IMMatrix[0]);
	}
	else
	{
		dKappa = atan(-IMMatrix[1] / IMMatrix[0]) + PI;
	}

	//计算镜头透视中心的位置,地心坐标转换为切面坐标系（即成图坐标系m）坐标
	THREEDPOINT curPoint;
	BLHToXYZ(dB, dL, dH, curPoint);
	double dXs = (curPoint.dX - XYZPoint.dX)*EMMatrix[0] + (curPoint.dY - XYZPoint.dY)*EMMatrix[1] + (curPoint.dZ - XYZPoint.dZ)*EMMatrix[2];
	double dYs = (curPoint.dX - XYZPoint.dX)*EMMatrix[3] + (curPoint.dY - XYZPoint.dY)*EMMatrix[4] + (curPoint.dZ - XYZPoint.dZ)*EMMatrix[5];
	double dZs = (curPoint.dX - XYZPoint.dX)*EMMatrix[6] + (curPoint.dY - XYZPoint.dY)*EMMatrix[7] + (curPoint.dZ - XYZPoint.dZ)*EMMatrix[8];



	fprintf(fEO, "\n%20.10f	%20.10f	%20.10f	%15.10f	%15.10f	%15.10f ", dXs, dYs, dZs, dPhi, dOmega, dKappa);
	for (int n = 0; n<9; n++)
	{
		fprintf(fEO, "	%10.4f", IMMatrix[n]);
	}

	return TRUE;
}




//=======================计算像元点坐标函数======================//
//      函数名称:  CalCoordinates
//		功能：计算所有影像地物像元（中间影像体）的所有大地坐标，并将计算得到的大地坐标保存到三维点数组中
//      参数：nHeight：		影像高度
//			  nImgWidth：	影像地物像元（中间影像体）
//			  nInterval：	影像EO文件行与影像行对应关系
//			  pGoundPt：	三维点数组
//      返回：BOOL
//      说明：
//===============================================================//
long  CalCoordinates(int nHeight, int nImgWidth, int nInterval, THREEDPOINT *pGoundPt, double *pdEOData, float fFov, float fFocalLen)
{
	int i = 0, j = 0;
	CPOINT temPt;

	double dPhotoPt[3];
	double dModelPt[3];
	memset(dPhotoPt, 0, 3 * sizeof(double));
	memset(dModelPt, 0, 3 * sizeof(double));


	double dAngle = 0;				//
	double dPixelIFov = double(fFov*PI / (180 * nImgWidth));
	unsigned long dwOffset = 0;
	double dXs = 0, dYs = 0, dZs = 0;

	double dRMatrix[9];
	for (i = 0; i<nHeight; i++)
	{
		dXs = pdEOData[i * 15];
		dYs = pdEOData[i * 15 + 1];
		dZs = pdEOData[i * 15 + 2];

		dRMatrix[0] = pdEOData[i * 15 + 6];
		dRMatrix[1] = pdEOData[i * 15 + 7];
		dRMatrix[2] = pdEOData[i * 15 + 8];
		dRMatrix[3] = pdEOData[i * 15 + 9];
		dRMatrix[4] = pdEOData[i * 15 + 10];
		dRMatrix[5] = pdEOData[i * 15 + 11];
		dRMatrix[6] = pdEOData[i * 15 + 12];
		dRMatrix[7] = pdEOData[i * 15 + 13];
		dRMatrix[8] = pdEOData[i * 15 + 14];
		for (j = 0; j<nImgWidth; j++)
		{
			temPt.x = -j + nImgWidth / 2;
			temPt.y = i;
			dAngle = temPt.x*dPixelIFov;

			memset(dPhotoPt, 0, 3 * sizeof(double));
			memset(dModelPt, 0, 3 * sizeof(double));
			dPhotoPt[0] = (fFocalLen)*sin(dAngle);
			dPhotoPt[1] = 0;
			dPhotoPt[2] = -(fFocalLen)*cos(dAngle);

			MatrixMuti(dRMatrix, dPhotoPt, dModelPt, 3, 3, 1);
			dwOffset = i*nImgWidth + j;
			pGoundPt[dwOffset].dX = dXs + (pGoundPt[dwOffset].dZ - dZs)*dModelPt[0] / dModelPt[2];
			pGoundPt[dwOffset].dY = dYs + (pGoundPt[dwOffset].dZ - dZs)*dModelPt[1] / dModelPt[2];
		}
	}

	return 0;
}
