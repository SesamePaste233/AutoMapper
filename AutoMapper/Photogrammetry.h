#pragma once
#include "CommonHeader.h"

namespace photo {

	class cameraParas {
	protected:
		cv::Mat camera_matrix;
		cv::Mat dist_coeffs;

	public:
		double factor;
		bool caliberated;

		cameraParas() {
			this->camera_matrix = cv::Mat::eye(3, 3, CV_64F);
			this->dist_coeffs = cv::Mat::zeros(1, 5, CV_64F);
			this->caliberated = false;
			this->factor = 1;
		};

		cameraParas(std::string FilePath) {
			readCameraParas(FilePath);
			caliberated = true;
		}

		cv::Mat cameraMatrix() {
			if (this->caliberated) {
				return this->camera_matrix;
			}
			return cv::Mat();
		};

		cv::Mat distCoeffs() {
			if (this->caliberated) {
				return this->dist_coeffs;
			}
			return cv::Mat();
		};

		double fx() {
			if (this->caliberated) {
				return this->camera_matrix.at<double>(0, 0);
			}
			return 0;
		}

		double fy() {
			if (this->caliberated) {
				return this->camera_matrix.at<double>(1, 1);
			}
			return 0;
		}

		double cx() {
			if (this->caliberated) {
				return this->camera_matrix.at<double>(0, 2);
			}
			return 0;
		}

		double cy() {
			if (this->caliberated) {
				return this->camera_matrix.at<double>(1, 2);
			}
			return 0;
		}

		void readCameraParas(std::string FilePath) {
			std::ifstream ifs(FilePath);
			float fx, fy, cx, cy, k1, k2, k3, p1, p2, fact;
			ifs >> fx >> fy >> cx >> cy >> k1 >> k2 >> p1 >> p2 >> k3 >> fact;
			this->camera_matrix = cv::Mat::eye(3, 3, CV_64F);
			this->dist_coeffs = cv::Mat::zeros(1, 5, CV_64F);
			this->camera_matrix.at<double>(0, 0) = fx / fact;
			this->camera_matrix.at<double>(1, 1) = fy / fact;
			this->camera_matrix.at<double>(0, 2) = cx / fact;
			this->camera_matrix.at<double>(1, 2) = cy / fact;
			this->dist_coeffs.at<double>(0) = k1;
			this->dist_coeffs.at<double>(1) = k2;
			this->dist_coeffs.at<double>(2) = p1;
			this->dist_coeffs.at<double>(3) = p2;
			this->dist_coeffs.at<double>(4) = k3;
			this->factor = fact;
			this->caliberated = true;
			ifs.close();
		}
	};

	template<int a, int b>
	class Pair {
	public:
		Eigen::RowVectorXd img_pt;
		Eigen::RowVectorXd obj_pt;
		Pair() {
			img_pt = Eigen::RowVectorXd(a);
			obj_pt = Eigen::RowVectorXd(b);
		}
		Pair(Eigen::RowVectorXd i, Eigen::RowVectorXd j) {
			img_pt = Eigen::RowVectorXd(a);
			obj_pt = Eigen::RowVectorXd(b);
			img_pt = i, obj_pt = j;
		}
	};

	Eigen::Matrix3d rotationMatrix(double r1, double r2, double r3);

	/*
	@brief This function calculates parameters of collinear equation based on the input.
	Warning: ***UNREALIZED***
	
	*/
	void collinearEquation(Eigen::VectorXd& intrinsic_elements, Eigen::VectorXd& extrinsicElems, Eigen::VectorXd& obj_pt, Eigen::VectorXd& img_pt);

	/*
	@brief Calculates affine transformation matrix. A * X' = X where X and X' is the homogeneous coordinates of image points.
	@param pairs: Storing points as X->img_pt, X'->obj_pt
	@param transformation_matrix:
	  A = 
	   {a1, a2, a0;
		b1, b2, b0;
		0,	0,	1;}
	@return x_error and y_error is the RMS of these points.
	*/
	Eigen::Vector2d interiorOrientation(std::vector<Pair<2, 2>> pairs, Eigen::Matrix3d& transformation_matrix);

	//@brief Calculates the extrinsic elements of a certain photo via paired image-object points
	Eigen::VectorXd calcExtrinsicElems(std::vector<Pair<2, 3>>& pairs, cameraParas& camera, double photographic_scale, opt::Method method = opt::Method::GaussNewton);

	//@brief Triangulates series of corresponding image points of the same object point.
	Eigen::Vector3d triangulatePoint(std::vector<Eigen::Vector2d>& img_pts, std::vector<Eigen::Vector3d>& extrinsic_linear_elems, std::vector<Eigen::Matrix3d>& rotationMatrixs, std::vector<cameraParas> & cameras);
	
	//@brief Calculates the 3D transformation between two cameras via paired image points.
	Eigen::VectorXd calcTransformation(std::vector<Pair<2, 2>>& pairs, Eigen::VectorXd& elements, cameraParas& camera1, cameraParas& camera2, opt::Method method = opt::Method::GaussNewton);
	
	//@brief Calculates the 3D transformation between model-coords and GCP-coords via paired 3D-points.
	Eigen::VectorXd calcTransformation(std::vector<Pair<3, 3>>& pairs, double& scale, Eigen::Matrix3d& rotation, Eigen::Vector3d& translation, bool simple_translation = true, opt::Method method = opt::Method::GaussNewton);

	Eigen::VectorXd sparseBundleAdjustment(cameraParas& camera, std::vector<std::map<int, Pair<2, 3>>>& pts, int num_img_pts, int num_obj_pts, std::vector<Eigen::VectorXd>& extrinsic_elements, opt::Method method = opt::Method::SparseLM);
}