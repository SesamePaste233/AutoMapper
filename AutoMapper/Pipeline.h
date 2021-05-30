#pragma once
/*


*/
#include "Photogrammetry.h"

namespace photo {
	class Pipeline {
	public:
		//Define storage structures for pipeline
		class Photo {
		public:
			std::string name;
			std::string camera_name;
			std::map<std::string, Eigen::RowVector2d> pts;
			bool interior_orientation_done;
			std::vector<Pair<2, 2>> inter_orient_pts;
			bool has_extrinsic_elems;
			Eigen::Vector3d extrinsicLinearElems;
			Eigen::Matrix3d rotationMatrix;

			//The project or zone THIS photo belongs to.
			Pipeline* master;

			//Next photo in this tie.
			Photo* next;

			Photo() :master(nullptr), next(nullptr) {
				interior_orientation_done = false;
				has_extrinsic_elems = false;
				camera_name = "default";
				extrinsicLinearElems.fill(0);
				rotationMatrix = Eigen::MatrixXd::Identity(3, 3);
			};

			Photo(Pipeline master) :master(&master), next(nullptr) {
				interior_orientation_done = false;
				has_extrinsic_elems = false;
				camera_name = "default";
				extrinsicLinearElems.fill(0);
				rotationMatrix = Eigen::MatrixXd::Identity(3, 3);
			};

			//Calculating extrinsic elements specific for THIS photo.
			bool calcExtrinsicElems(bool verbose = false);

			bool interiorOrientate(bool verbose = false);

			void applyTransformation(Eigen::VectorXd trans_vec, double scale = 1);

			void applyTransformation(Eigen::Vector3d tvec,Eigen::Matrix3d rotation, double scale = 1);
		};

		class PhotogrammetricPoint {
		public:
			std::string name;
			std::vector<std::string> in_photos;
			std::map<int, std::vector<std::string>> in_tie;
			Eigen::Vector3d coordinates;

			PhotogrammetricPoint() {};
			PhotogrammetricPoint(std::string name):name(name) {
				coordinates.fill(0);
			};
		};

		class TieModel {
		public:
			int tie_id;
			bool is_free;
			std::vector<PhotogrammetricPoint> model_pts;
			TieModel():is_free(true) {};
		};

		class zoneInfo {
		public:
			double photographic_scale;

			std::vector<std::vector<std::string>>ties;

			zoneInfo() :photographic_scale(10000) {
			};
			zoneInfo(std::string zone_info_file) :photographic_scale(10000) {
				this->readZoneInfo(zone_info_file);
			};
			double m() {
				return photographic_scale;
			}
			void readZoneInfo(std::string zone_info_file) {
				std::ifstream ifs(zone_info_file);
				while (!ifs.eof()) {
					double m;
					ifs >> m;
					this->photographic_scale = m;
					int name;
					std::vector<std::string> tie;
					while (!ifs.eof()) {
						ifs >> name;
						if (std::to_string(name) == "-99") {
							break;
						}
						tie.push_back(std::to_string(name));
					}
					if (!tie.empty()) {
						ties.push_back(tie);
					}
				}
			}
		};

		//Information of the object zone
		zoneInfo zone;

		//Tables and tuples for storing cameras, photos and points. Name-value pair.
		std::map<std::string, cameraParas> cameras;

		std::map<std::string, Photo> photos;

		std::map<std::string, Eigen::RowVector3d> gcps;

		std::map<std::string, PhotogrammetricPoint> obj_pts;

		std::vector<TieModel> tie_models;

		//Status indicators
		bool zone_info_set = false;

		bool camera_info_set = false;

		Pipeline() {};
		Pipeline(std::string zone_info_file, std::string img_pts_file, std::string gcp_file) {
			addZoneInfo(zone_info_file);
			readImagePoints(img_pts_file);
			readControlPoints(gcp_file);
		};

		void addZoneInfo(std::string zone_info_file) {
			zone = zoneInfo(zone_info_file);
			zone_info_set = true;
			arrangeTies();
		}

		void addZoneInfo(zoneInfo zone_info) {
			this->zone = zone_info;
			zone_info_set = true;
			arrangeTies();
		}

		void addCameraInfo(std::string name, cameraParas camera) {
			cameras[name] = camera;
			camera_info_set = true;
		}

		void addCameraInfo(std::string name, std::string camera_file) {
			cameras[name] = cameraParas(camera_file);
			camera_info_set = true;
		}

		void readImagePoints(std::string img_pts_file);

		void readInterOrientFile(std::string inter_orient_file);

		void readControlPoints(std::string gcp_file);

		void applyInteriorOrientation(bool verbose = false);

		void calcExtrinsicElemsByPhotos(bool verbose = false);

		//calculate transformation from this photo to the next photo
		void calcTransformationByPhotos(std::string photo_name, bool verbose = false);

		void calcTransformationByTie(int tie_index, bool verbose = false);

		void triangulateAllImagePoints(bool triangulate_gcps = true, bool verbose = false);

		void triangulateImagePointsByTie(int tie_id, bool triangulate_gcps = true, bool verbose = false);

		void formZoneNet();

		inline void printSeparator() {
			std::cout << "===========================" << std::endl;
		};

	protected:

		void arrangeTies();

		void mergeModels();
	};
	
};