#include "opencv2/opencv.hpp"
#include "opencl.h"

int main(int argc, char **argv)
{

	cv::VideoCapture cap(0);
	cv::Mat frame;
	cv::Mat result(480, 640, CV_8UC1);

	OpenCL::initialize(0, 0);

	if (!cap.isOpened())
		return -1;

	while (1)
	{
		do
		{
			cap >> frame;
		} while (frame.empty());

		OpenCL::detectLine(result.data, frame.data);

		cv::imshow("test", result);

		int key = cv::waitKey(1);
		if (key == 113)
			break;

	}

	OpenCL::release();

}