#pragma once
#include "CommonHeader.h"

//Used types. Mainly written out of security issues.
namespace types {
	class Array {
	public:
		typedef const cv::Mat& Input;
		typedef cv::Mat& Output, InputOutput;
	};

	typedef Array::Input& InputArray;
	typedef Array::Output& OutputArray;
	typedef Array::InputOutput& InputOutputArray;


	template<class T>
	class Sequence {
		Sequence() = delete;
	public:
		typedef const std::vector<T>& Input;
		typedef std::vector<T>& Output, InputOutput;
	};

	typedef types::Sequence<cv::Point2f>::Input& InputSeqOfPoints;
	typedef types::Sequence<cv::Point2f>::Output& OutputSeqOfPoints;
	typedef types::Sequence<cv::Point2f>::InputOutput& InputOutputSeqOfPoints;

	typedef types::Sequence<cv::Point3f>::Input& InputSeqOf3DPoints;
	typedef types::Sequence<cv::Point3f>::Output& OutputSeqOf3DPoints;
	typedef types::Sequence<cv::Point3f>::InputOutput& InputOutputSeqOf3DPoints;


	//Eqivalent to std::vector<cv::Mat>
	typedef types::Sequence<cv::Mat>::Input& InputSeqOfArrays;
	typedef types::Sequence<cv::Mat>::Output& OutputSeqOfArrays;
	typedef types::Sequence<cv::Mat>::InputOutput& InputOutputSeqOfArrays;

}