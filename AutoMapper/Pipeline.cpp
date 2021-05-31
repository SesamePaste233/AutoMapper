#include "Pipeline.h"

void photo::Pipeline::readImagePoints(std::string img_pts_file)
{
    std::ifstream ifs(img_pts_file);
    while (!ifs.eof()) {
        int photo_name, camera_id;
        ifs >> photo_name >> camera_id;
        photos[std::to_string(photo_name)].name = std::to_string(photo_name);
        photos[std::to_string(photo_name)].camera_name = std::to_string(camera_id);
        photos[std::to_string(photo_name)].master = this;
        int name;
        double x, y;
        while (!ifs.eof()) {
            ifs >> name >> x >> y;
            if (std::to_string(name) == "-99") {
                break;
            }
            photos[std::to_string(photo_name)].pts[std::to_string(name)] = Eigen::RowVector2d(x, y);
            obj_pts[std::to_string(name)].name=std::to_string(name);
            obj_pts[std::to_string(name)].in_photos.push_back(std::to_string(photo_name));
            for (int i = 0;i < zone.ties.size();i++) {
                for (auto phot : zone.ties[i]) {
                    if (std::to_string(photo_name) == phot) {
                        obj_pts[std::to_string(name)].in_tie[i].push_back(phot);
                    }
                }
            }
        }
    }

}

void photo::Pipeline::readInterOrientFile(std::string inter_orient_file)
{
    std::ifstream ifs(inter_orient_file);
    while (!ifs.eof()) {
        int photo_name;
        ifs >> photo_name;
        photos[std::to_string(photo_name)].name = std::to_string(photo_name);
        photos[std::to_string(photo_name)].master = this;
        int id;
        double x, y, _x, _y;
        while (!ifs.eof()) {
            ifs >> id;
            if (id == -99) {
                break;
            } 
            ifs >> x >> y >> _x >> _y;
            Pair<2, 2> pair;
            pair.obj_pt << x, y;
            pair.img_pt << _x, _y;
            photos[std::to_string(photo_name)].inter_orient_pts.push_back(pair);
        }
    }
}

void photo::Pipeline::readControlPoints(std::string gcp_file)
{
    std::ifstream ifs(gcp_file);
    while (!ifs.eof()) {
        int name;
        double x, y, z;
        ifs >> name >> x >> y >> z;
        gcps[std::to_string(name)] = Eigen::RowVector3d(x, y, z);
    }
}

void photo::Pipeline::applyInteriorOrientation(bool verbose)
{
    printSeparator();
    std::cout << "INFO: Applying interior orientation to each photo." << std::endl;
    int n = 0;
    for (auto iter = photos.begin();iter != photos.end();iter++) {
        std::cout << "INFO: Calculating photo: " + iter->first<< "(" << ++n << "/" << photos.size() << ") ..." << std::endl;
        clock_t t1 = clock();
        if (!iter->second.interiorOrientate(verbose)) {
            std::cout << "INFO: Failed to orient photo: " + iter->first + ". Skipping... " << std::endl;
            continue;
        }
        clock_t t2 = clock();
        std::cout << "INFO: Calulation finished for photo: " + iter->first + ". Elapsed time: " + std::to_string((t2 - t1) / float(CLOCKS_PER_SEC) * 1000) + " ms. " << std::endl;
    }
}

void photo::Pipeline::calcExtrinsicElemsByPhotos(bool verbose)
{
    printSeparator();
    if (!zone_info_set || !camera_info_set) {
        std::cout << "WARN: Please set zone informations and add cameras before solving." << std::endl;
        return;
    }
    std::cout << "INFO: Calculating extrinsic elements of each photo." << std::endl;
    int n = 0;
    for (auto iter = photos.begin();iter != photos.end();iter++) {
        std::cout << "INFO: Calculating photo: " + iter->first << "(" << ++n << "/" << photos.size() << ") ..." << std::endl;
        clock_t t1 = clock();
        if (!iter->second.calcExtrinsicElems(verbose)) {
            std::cout << "INFO: No enough gcps for photo: " + iter->first + ". Skipping... " << std::endl;
            continue;
        }
        clock_t t2 = clock();
        std::cout << "INFO: Calulation finished for photo: " + iter->first + ". Elapsed time: " + std::to_string((t2 - t1) / float(CLOCKS_PER_SEC) * 1000) + " ms. " << std::endl;
    }
}

void photo::Pipeline::calcTransformationByPhotos(std::string photo_name, bool verbose)
{
    printSeparator();
    if (photos[photo_name].next == nullptr) {
        std::cout << "INFO: Reaching the end of the tie. Finished." << std::endl;
        return;
    }
    std::cout << "INFO: Calculating transformation from photo: "<<photo_name<<" to "<<photos[photo_name].next->name<<" ." << std::endl;

    clock_t t1 = clock();

    std::vector<Pair<2, 2>> pairs;
    Eigen::VectorXd elements(5);
    cameraParas camera1,camera2;
    try {
        auto c_iter = photos[photo_name].master->cameras.find(photos[photo_name].camera_name);
        if (c_iter != photos[photo_name].master->cameras.end()) {
            camera1 = c_iter->second;
        }
        else {
            camera1 = photos[photo_name].master->cameras.begin()->second;
        }
        c_iter = photos[photo_name].next->master->cameras.find(photos[photo_name].next->camera_name);
        if (c_iter != photos[photo_name].next->master->cameras.end()) {
            camera2 = c_iter->second;
        }
        else {
            camera2 = photos[photo_name].next->master->cameras.begin()->second;
        }
    }
    catch (...) {
        throw("ERROR: At photo::Pipeline::calcTransformationByPhotos(). Bad assertion. Please check if these photos has the right master. ");
    }
    for (auto pt1 = photos[photo_name].pts.begin();pt1 != photos[photo_name].pts.end();pt1++) {
        for (auto pt2 = photos[photo_name].next->pts.begin();pt2 != photos[photo_name].next->pts.end();pt2++) {
            if (pt1->first == pt2->first) {
                Pair<2, 2> pair;
                pair.img_pt = pt1->second;
                pair.obj_pt = pt2->second;
                pairs.push_back(pair);
                break;
            }
        }
    }

    auto error = photo::calcTransformation(pairs, elements, camera1, camera2);

    Eigen::Vector3d translation;
    Eigen::Matrix3d rotation;
    translation << 1, tan(elements(3)), (1 + tan(elements(3))) * tan(elements(4));
    translation.normalize();
    rotation = photo::rotationMatrix(elements(0), elements(1), elements(2));
    photos[photo_name].next->applyTransformation(translation, rotation, 1);

    if (verbose) {
        std::cout << " Translation (normalized) : " << std::endl << translation << std::endl;
        std::cout << " Rotation matrix: " << std::endl << rotation << std::endl;
        std::cout << " Error (norm2 q): " << sqrt(error.norm()*error.norm()/error.rows()) * 1000 << "mm. " << std::endl;
    }

    clock_t t2 = clock();
    std::cout << "INFO: Calulation finished. Elapsed time: " + std::to_string((t2 - t1) / float(CLOCKS_PER_SEC) * 1000) + " ms. " << std::endl;
}

void photo::Pipeline::calcTransformationByTie(int tie_index, bool verbose)
{
    printSeparator();
    if (tie_index < 0 || tie_index >= zone.ties.size()) {
        std::cout << "ERROR: Index out of range." << std::endl;
        return;
    }
    std::cout << "INFO: Calculating transformating with tie id: " << tie_index << "." << std::endl;
    for (int i = 0;i < zone.ties[tie_index].size() - 1;i++) {
        this->calcTransformationByPhotos(zone.ties[tie_index][i], verbose);
        if (i > 0) {
            std::vector<std::vector<Eigen::Vector2d>> pair12, pair23;
            std::vector<Eigen::Vector3d> le12, le23;
            std::vector<Eigen::Matrix3d> r12, r23;
            std::vector<cameraParas> c12, c23;
            for (auto pt1 : photos[zone.ties[tie_index][i - 1]].pts) {
                for (auto pt2 : photos[zone.ties[tie_index][i]].pts) {
                    if (pt1.first == pt2.first) {
                        for (auto pt3 : photos[zone.ties[tie_index][i + 1]].pts) {
                            if (pt2.first == pt3.first) {
                                std::vector<Eigen::Vector2d> pt_pair12, pt_pair23;
                                pt_pair12.push_back(pt1.second), pt_pair12.push_back(pt2.second);
                                pair12.push_back(pt_pair12);
                                pt_pair23.push_back(pt2.second), pt_pair23.push_back(pt3.second);
                                pair23.push_back(pt_pair23);
                            }
                        }
                    }
                }
            }
            auto c_iter = photos[zone.ties[tie_index][i - 1]].master->cameras.find(photos[zone.ties[tie_index][i - 1]].camera_name);
            if (c_iter != photos[zone.ties[tie_index][i - 1]].master->cameras.end()) {
                c12.push_back(c_iter->second);
            }
            else {
                c12.push_back(photos[zone.ties[tie_index][i - 1]].master->cameras.begin()->second);
            }
            c_iter = photos[zone.ties[tie_index][i]].next->master->cameras.find(photos[zone.ties[tie_index][i]].next->camera_name);
            if (c_iter != photos[zone.ties[tie_index][i]].next->master->cameras.end()) {
                c12.push_back(c_iter->second);
                c23.push_back(c_iter->second);
            }
            else {
                c12.push_back(photos[zone.ties[tie_index][i]].next->master->cameras.begin()->second);
                c23.push_back(photos[zone.ties[tie_index][i]].next->master->cameras.begin()->second);
            }
            c_iter = photos[zone.ties[tie_index][i + 1]].next->master->cameras.find(photos[zone.ties[tie_index][i + 1]].next->camera_name);
            if (c_iter != photos[zone.ties[tie_index][i + 1]].next->master->cameras.end()) {
                c23.push_back(c_iter->second);
            }
            else {
                c23.push_back(photos[zone.ties[tie_index][i + 1]].next->master->cameras.begin()->second);
            }

            le12.push_back(photos[zone.ties[tie_index][i - 1]].extrinsicLinearElems);
            le12.push_back(photos[zone.ties[tie_index][i]].extrinsicLinearElems);
            le23.push_back(photos[zone.ties[tie_index][i]].extrinsicLinearElems);
            le23.push_back(photos[zone.ties[tie_index][i + 1]].extrinsicLinearElems);
            r12.push_back(photos[zone.ties[tie_index][i - 1]].rotationMatrix);
            r12.push_back(photos[zone.ties[tie_index][i]].rotationMatrix);
            r23.push_back(photos[zone.ties[tie_index][i]].rotationMatrix);
            r23.push_back(photos[zone.ties[tie_index][i + 1]].rotationMatrix);

            std::vector<Eigen::Vector3d> obj_pts1, obj_pts2;
            for (int j = 0;j < pair12.size();j++) {
                obj_pts1.push_back(triangulatePoint(pair12[i], le12, r12, c12));
                obj_pts2.push_back(triangulatePoint(pair23[i], le23, r23, c23));
            }
            double scale = 0;
            for (int j = 0;j < obj_pts1.size();j++) {
                scale += (obj_pts1[j](2) - photos[zone.ties[tie_index][i]].extrinsicLinearElems(3)) / obj_pts2[j](2);
            }
            photos[zone.ties[tie_index][i + 1]].applyTransformation(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity(), scale);
        }
    }
}

void photo::Pipeline::triangulateAllImagePoints(bool triangulate_gcps, bool verbose)
{
    printSeparator();
    for (auto iter = obj_pts.begin();iter != obj_pts.end();iter++) {
        if (iter->second.in_photos.size() > 1) {//At least two corresponding image-points are needed to solve
            if (!triangulate_gcps && this->gcps.find(iter->first) != gcps.end()) {//if this point is a gcp
                continue;
            }
            std::vector<Eigen::Vector2d> pts;
            std::vector<Eigen::Vector3d> extern_paras;
            std::vector<Eigen::Matrix3d> rotationMatrixs;
            std::vector<cameraParas> cameras;
            //search for photo containing this image point
            for (int i = 0;i < iter->second.in_photos.size();i++) {
                //constructing parameters for solving

                if (!photos[iter->second.in_photos[i]].has_extrinsic_elems) {
                    //if this photo has none calculated extrinsic elements. Skip this one.
                    std::cout << "WARN: Image point: " << iter->first << " in photo: " << iter->second.in_photos[i] << " which has no trusted extrinsic elemsnts. Calculate the extrinsic elements for this photo before taking this image point into triangulation." << std::endl;
                    continue;
                }

                pts.push_back(photos[iter->second.in_photos[i]].pts[iter->first].transpose());

                extern_paras.push_back(photos[iter->second.in_photos[i]].extrinsicLinearElems);

                rotationMatrixs.push_back(photos[iter->second.in_photos[i]].rotationMatrix);

                cameraParas camera;
                auto c_iter = this->cameras.find(photos[iter->second.in_photos[i]].camera_name);
                if (c_iter != this->cameras.end()) {
                    camera = c_iter->second;
                }
                else {
                    camera = this->cameras.begin()->second;
                }
                cameras.push_back(camera);
            }
            //solve the coords of this point
            iter->second.coordinates = photo::triangulatePoint(pts, extern_paras, rotationMatrixs, cameras);
        }
        else {
            std::cout << "WARN: Image point: " << iter->first << " appears to be only in " << iter->second.in_photos.size() << " photo, which is insufficient for calculation. Skipping..." << std::endl;
        }
    }
}

void photo::Pipeline::triangulateImagePointsByTie(int tie_id, bool triangulate_gcps, bool verbose)
{
    for (auto obj_pt : obj_pts) {
        auto _photos = obj_pt.second.in_tie.find(tie_id);
        if (_photos != obj_pt.second.in_tie.end()) {
            if (_photos->second.size() < 2) {
                continue;
            }
            std::vector<Eigen::Vector2d> pts;
            std::vector<Eigen::Vector3d> extern_paras;
            std::vector<Eigen::Matrix3d> rotationMatrixs;
            std::vector<cameraParas> cameras;
            //search for photo containing this image point
            for (int i = 0;i < _photos->second.size();i++) {
                //constructing parameters for solving
                pts.push_back(photos[_photos->second[i]].pts[obj_pt.first].transpose());

                extern_paras.push_back(photos[_photos->second[i]].extrinsicLinearElems);

                rotationMatrixs.push_back(photos[_photos->second[i]].rotationMatrix);

                cameraParas camera;
                auto c_iter = this->cameras.find(photos[_photos->second[i]].camera_name);
                if (c_iter != this->cameras.end()) {
                    camera = c_iter->second;
                }
                else {
                    camera = this->cameras.begin()->second;
                }
                cameras.push_back(camera);
            }
            //solve the coords of this point
            obj_pt.second.coordinates = photo::triangulatePoint(pts, extern_paras, rotationMatrixs, cameras);
            tie_models[tie_id].model_pts.push_back(obj_pt.second);
        }
    }
}

void photo::Pipeline::formZoneNet()
{
    for (int i = 0;i < zone.ties.size();i++) {
        calcTransformationByTie(i);
        triangulateImagePointsByTie(i);
    }
    mergeModels();
    triangulateAllImagePoints();
}


void photo::Pipeline::arrangeTies()
{
    for (auto photo:photos) {
        photo.second.next = nullptr;
    }
    for (auto tie : zone.ties) {
        for (auto iter2 = tie.begin();iter2 + 1 != tie.end();iter2++) {
            this->photos[*iter2].master = this;
            this->photos[*(iter2 + 1)].master = this;
            this->photos[*iter2].next = &this->photos[*(iter2 + 1)];
        }
    }
}

void photo::Pipeline::mergeModels()
{
    for (int i = 0;i < tie_models.size()-1;i++) {
        std::vector<Pair<3, 3>>pairs;
        for (auto pt1 : tie_models[i].model_pts) {
            for (auto pt2 : tie_models[i + 1].model_pts) {
                if (pt1.name == pt2.name) {
                    pairs.push_back(Pair<3, 3>(pt1.coordinates, pt2.coordinates));
                }
            }
        }
        if (pairs.size() < 3) {
            throw;
        }
        double scale=0;
        Eigen::Vector3d translation;
        Eigen::Matrix3d rotation;
        translation.fill(0);
        rotation = Eigen::Matrix3d::Identity();
        photo::calcTransformation(pairs, scale, rotation, translation);
        photos[zone.ties[i + 1][0]].applyTransformation(translation, rotation, scale);
    }
}

bool photo::Pipeline::Photo::calcExtrinsicElems(bool verbose)
{
    //find corresponding camera
    cameraParas camera;
    try {
        auto c_iter = master->cameras.find(camera_name);
        if (c_iter != master->cameras.end()) {
            camera = c_iter->second;
        }
        else {
            camera = master->cameras.begin()->second;
        }
    }catch (...) {
        throw("ERROR: At photo::Pipeline::Photo::calcExtrinsicElems(). Bad assertion. Please check if THIS photo has the right master. ");
    }
    //find and organize image points and their corresponding object points
    std::vector<Pair<2, 3>> pairs;
    auto p_iter = pairs.begin();
    for (auto iter1 = pts.begin();iter1 != pts.end();iter1++) {
        auto iter2 = master->gcps.find(iter1->first);
        if (iter2 != master->gcps.end()) {
            Pair<2, 3> p;
            p.img_pt = iter1->second/1000, p.obj_pt = iter2->second;
            pairs.push_back(p);
        }
    }
    if (pairs.size() < 3) {
        return false;
    }
    //begin calculation
    auto extrinsicElems = photo::calcExtrinsicElems(pairs, camera, master->zone.m());
    this->extrinsicLinearElems = extrinsicElems.head(3);
    this->rotationMatrix = photo::rotationMatrix(extrinsicElems(3), extrinsicElems(4), extrinsicElems(5));
    this->has_extrinsic_elems = true;

    if (verbose) {
        std::cout << " Linear elements:" << std::endl;
        std::cout << extrinsicLinearElems << std::endl;
        std::cout << " Rotation matrix:" << std::endl;
        std::cout << rotationMatrix << std::endl;
    }
    return true;
}

bool photo::Pipeline::Photo::interiorOrientate(bool verbose)
{
    if (this->interior_orientation_done == true) {
        return true;
    }
    Eigen::Matrix3d trans;
    auto errors = photo::interiorOrientation(this->inter_orient_pts, trans);
    if (errors(0) == 0 && errors(1) == 0) {
        return false;
    }
    else {
        if (verbose) { 
            std::cout << " X error: " << errors(0) << ". " << std::endl << " Y error: " << errors(1) << ". " << std::endl;
        }
        //Apply transformation matrix to all image points in this photo.
        for (auto iter = this->pts.begin();iter != pts.end();iter++) {
            Eigen::Vector3d homo_pt;
            homo_pt.fill(1);
            homo_pt.head(2) = iter->second.transpose();
            iter->second = (trans.inverse() * homo_pt).head(2).transpose()/1000;
        }
        this->interior_orientation_done = true;
    }
    return true;
}

void photo::Pipeline::Photo::applyTransformation(Eigen::VectorXd trans_vec, double scale)
{
    this->rotationMatrix = photo::rotationMatrix(trans_vec(3), trans_vec(4), trans_vec(5))*this->rotationMatrix;
    this->extrinsicLinearElems += trans_vec.head(3);
    if (next == nullptr) {
        return;
    }
    else {
        this->next->applyTransformation(trans_vec, scale);
    }
}

void photo::Pipeline::Photo::applyTransformation(Eigen::Vector3d tvec, Eigen::Matrix3d rotation, double scale)
{
    this->rotationMatrix = rotation * this->rotationMatrix;
    this->extrinsicLinearElems += tvec;
    if (next == nullptr) {
        return;
    }
    else {
        this->next->applyTransformation(tvec, rotation, scale);
    }
}
