#include "Photogrammetry.h"


void photo::collinearEquation(Eigen::VectorXd& intrinsic_elements, Eigen::VectorXd& extrinsicElems, Eigen::VectorXd& obj_pt, Eigen::VectorXd& img_pt)
{
    
}

Eigen::Vector2d photo::interiorOrientation(std::vector<Pair<2, 2>> pairs, Eigen::Matrix3d& transformation_matrix)
{
    transformation_matrix = Eigen::MatrixXd::Identity(3, 3);
    if (pairs.size() <= 3) {
        std::cout << "ERROR: No enough points for calculation." << std::endl;
        return Eigen::Vector2d(0,0);
    }
    Eigen::MatrixXd A(pairs.size(), 3);
    Eigen::VectorXd Lx(pairs.size());
    Eigen::VectorXd Ly(pairs.size());
    A.fill(1);
    for (int i = 0;i < pairs.size();i++) {
        A.row(i).tail(2) << pairs[i].obj_pt(0), pairs[i].obj_pt(1);
        Lx(i) = pairs[i].img_pt(0);
        Ly(i) = pairs[i].img_pt(1);
    }

    auto a = (A.transpose() * A).inverse() * A.transpose() * Lx;
    auto b = (A.transpose() * A).inverse() * A.transpose() * Ly;
    transformation_matrix << a(1), a(2), a(0),
        b(1), b(2), b(0),
        0, 0, 1;

    Eigen::VectorXd x_error = (A * a - Lx);
    Eigen::VectorXd y_error = (A * b - Ly);

    return Eigen::Vector2d(sqrt(x_error.norm() * x_error.norm() / x_error.rows()),sqrt(y_error.norm() * y_error.norm() / y_error.rows()));
}

Eigen::VectorXd photo::calcExtrinsicElems(std::vector<Pair<2, 3>>& pairs, cameraParas& camera, double photographic_scale, opt::Method method)
{
    int x0 = camera.cx();
    int y0 = camera.cy();
    double f = camera.fx();
    Eigen::VectorXd X(6);//extrinsic elements
    X.fill(0);
    Eigen::MatrixXd A(2 * pairs.size(), 6);//Jacobian matrix
    Eigen::VectorXd L(2 * pairs.size());
    A.fill(0);
    L.fill(0);

    for (auto pair : pairs) {
        X(0) += pair.obj_pt[0];
        X(1) += pair.obj_pt[1];
    }
    X(0) /= pairs.size();
    X(1) /= pairs.size();
    X(2) = photographic_scale * f;

    opt::LeastSquareSolver solver(A, L, X);//declare a least-square solver as optimizer

    //constructing object function for optimization
    auto obj_fcn = [&](void*)->bool {
        for (int n = 0;n < pairs.size();n++) {
            Eigen::Matrix3d R = rotationMatrix(X(3), X(4), X(5));
            double Xslash = (pairs[n].obj_pt - X.head(3).transpose()) * R.col(0);
            double Yslash = (pairs[n].obj_pt - X.head(3).transpose()) * R.col(1);
            double Zslash = (pairs[n].obj_pt - X.head(3).transpose()) * R.col(2);
            double x = pairs[n].img_pt[0];
            double y = pairs[n].img_pt[1];
            L(0 + n * 2) = x - x0 + f * Xslash / Zslash;
            L(1 + n * 2) = y - y0 + f * Yslash / Zslash;
            for (int i = 0;i < 3;i++) {
                A(0 + n * 2, i) = (R(i, 0) * f + R(i, 2) * (x - x0)) / Zslash;
                A(1 + n * 2, i) = (R(i, 1) * f + R(i, 2) * (y - y0)) / Zslash;
            }
            A(0 + n * 2, 3) = (y - y0) * sin(X(4)) - ((x - x0) / f * ((x - x0) * cos(X(5)) - (y - y0) * sin(X(5))) + f * cos(X(5))) * cos(X(4));
            A(1 + n * 2, 3) = -(x - x0) * sin(X(4)) - ((y - y0) / f * ((x - x0) * cos(X(5)) - (y - y0) * sin(X(5))) - f * sin(X(5))) * cos(X(4));
            A(0 + n * 2, 4) = -f * sin(X(5)) - (x - x0) / f * ((x - x0) * sin(X(5)) + (y - y0) * cos(X(5)));
            A(1 + n * 2, 4) = -f * cos(X(5)) - (y - y0) / f * ((x - x0) * sin(X(5)) + (y - y0) * cos(X(5)));
            A(0 + n * 2, 5) = y - y0;
            A(1 + n * 2, 5) = -(x - x0);
        }
        return true;
    };

    solver.setObjectFunction(obj_fcn);

    solver.maxIteration() = 40;

    solver.method = method;

    auto flag = solver.solve(nullptr);
    
    if (flag == opt::QuitFlag::MaxIterReached) {
        std::cout << "ERROR: Max iteration reached. Result unfeasable." << std::endl;
        return Eigen::VectorXd(6);
    }

    return X;
}

Eigen::Vector3d photo::triangulatePoint(std::vector<Eigen::Vector2d>& img_pts, std::vector<Eigen::Vector3d>& extrinsic_linear_elems, std::vector<Eigen::Matrix3d>& rotationMatrixs, std::vector<cameraParas> & cameras)
{
    if (img_pts.size() <= 1)return Eigen::Vector3d();
    Eigen::MatrixX3d A(2 * img_pts.size(),3);
    Eigen::VectorXd L(2 * img_pts.size());
    for (int n = 0;n < img_pts.size();n++) {
        auto R = rotationMatrixs[n];
        double x = img_pts[n][0];
        double y = img_pts[n][1];
        double f = cameras[n].fx();
        double x0 = cameras[n].cx();
        double y0 = cameras[n].cy();
        for (int i = 0;i < 3;i++) {
            A(0 + n * 2, i) = R(i, 0) * f + R(i, 2) * (x - x0);
            A(1 + n * 2, i) = R(i, 1) * f + R(i, 2) * (y - y0);
        }
        for (int i = 0;i < 3;i++) {
            L(0 + n * 2) += extrinsic_linear_elems[n][i] * A(0 + n * 2, i);
            L(1 + n * 2) += extrinsic_linear_elems[n][i] * A(1 + n * 2, i);
        }
    }
    return (A.transpose() * A).inverse() * A.transpose() * L;
}

Eigen::VectorXd photo::calcTransformation(std::vector<Pair<2, 2>>& pairs, Eigen::VectorXd& elements, cameraParas & camera1, cameraParas& camera2, opt::Method method)
{
    int x0_1 = camera1.cx();
    int y0_1 = camera1.cy();
    double f_1 = camera1.fx();
    int x0_2 = camera2.cx();
    int y0_2 = camera2.cy();
    double f_2 = camera2.fx();
    double Bx = 1;

    Eigen::VectorXd X(5);//extrinsic elements
    Eigen::MatrixXd A(pairs.size(), 5);//Jacobian matrix
    Eigen::VectorXd L(pairs.size());
    X.fill(0);

    Eigen::MatrixXd pts1(pairs.size(), 3);
    for (int i = 0;i < pairs.size();i++) {
        pts1.row(i) << pairs[i].img_pt(0), pairs[i].img_pt(1), -f_1;
    }

    opt::LeastSquareSolver solver(A, L, X);

    auto obj_fcn = [&](void*)->bool {
        auto R_L = rotationMatrix(type, X(0), X(1), X(2));
        double By = Bx * tan(X(3));
        double Bz = sqrt(Bx * Bx + By * By) * tan(X(4));
        for (int j = 0;j < pairs.size();j++) {
            Eigen::Vector3d pts2;
            pts2 << pairs[j].obj_pt(0), pairs[j].obj_pt(1), -f_2;
            pts2 = R_L * pts2;
            double N1 = (Bx * pts2(2) - Bz * pts2(0)) / (pts1(j, 0) * pts2(2) - pts1(j, 2) * pts2(0));
            double N2 = (Bx * pts1(j, 2) - Bz * pts1(j, 0)) / (pts1(j, 0) * pts2(2) - pts1(j, 2) * pts2(0));
            A(j, 0) = -pts2(0) * pts2(1) / pts2(2) * N2;
            A(j, 1) = -(pts2(2) + pts2(1) * pts2(1) / pts2(2)) * N2;
            A(j, 2) = pts2(0) * N2;
            A(j, 3) = Bx;
            A(j, 4) = -pts2(1) / pts2(2) * Bx;
            L(j) = N1 * pts1(j, 1) - N2 * pts2(1) - By;
        }
        return true;
    };

    solver.setObjectFunction(obj_fcn);

    solver.method = method;

    auto flag = solver.solve(nullptr);

    if (flag == opt::QuitFlag::MaxIterReached) {
        std::cout << "ERROR: Max iteration reached. Result unfeasable." << std::endl;
        return Eigen::VectorXd(5);
    }

    elements = X;

    while (!obj_fcn(nullptr));

    return A*X-L;
}

Eigen::VectorXd photo::calcTransformation(std::vector<Pair<3, 3>>& pairs, double& scale, Eigen::Matrix3d& rotation, Eigen::Vector3d& translation, bool simple_translation, opt::Method method)
{
    if (pairs.size() < 3) {
        return Eigen::VectorXd();
    }
    if (simple_translation) {
        Eigen::VectorXd X(4);//lambda, phi, omega, kappa
        Eigen::MatrixXd A(3, 4);
        Eigen::VectorXd L(3);
        X.fill(0);
        X(0) = 1;
        Eigen::Vector3d x1_mass, x2_mass;
        x1_mass.fill(0);
        x2_mass.fill(0);
        for (auto pair : pairs) {
            x1_mass += pair.img_pt;
            x2_mass += pair.obj_pt;
        }
        x1_mass /= pairs.size();
        x2_mass /= pairs.size();

        int pt_count = 0;

        auto obj_fcn = [&](void*)->bool {
            Eigen::Vector3d Xhat, _Xhat;
            Xhat = pairs[pt_count].img_pt.transpose() - x1_mass;
            _Xhat = pairs[pt_count].obj_pt.transpose() - x2_mass;
            A.col(0) = Xhat.transpose();
            A.col(1) << -X(0) * Xhat(2), 0, X(0)* Xhat(0);
            A.col(2) << -Xhat(1) * sin(X(1)), Xhat(0)* sin(X(1)) - Xhat(2) * cos(X(1)), Xhat(1)* cos(X(1));
            A.col(2) *= X(0);
            A.col(3) << -Xhat(1) * cos(X(1)) * cos(X(2)) - Xhat(2) * sin(X(2)),
                Xhat(0)* cos(X(1))* cos(X(2)) + Xhat(2) * sin(X(1)) * cos(X(2)),
                Xhat(0)* sin(X(2)) - Xhat(1) * sin(X(1)) * cos(X(2));
            A.col(3) *= X(0);
            L = _Xhat - rotationMatrix(X(1), X(2), X(3)) * Xhat * X(0);

            if (++pt_count >= pairs.size()) {
                pt_count = 0;
                return true;
            }
            else {
                return false;
            }
        };

        opt::LeastSquareSolver solver(A, L, X);

        solver.setObjectFunction(obj_fcn);

        solver.method = method;

        solver.function_num = 3 * pairs.size();

        auto flag = solver.solve(nullptr);

        if (flag == opt::QuitFlag::MaxIterReached) {
            std::cout << "ERROR: Max iteration reached. Result unfeasable." << std::endl;
            return Eigen::VectorXd(3);
        }

        rotation = rotationMatrix(X(1), X(2), X(3));
        scale = X(0);
        translation = x2_mass - X(0) * rotationMatrix(PhyOmegaKappa, X(1), X(2), X(3)) * x1_mass;
    }
    else {
        Eigen::VectorXd X(7);//lambda, phi, omega, kappa,x,y,z
        Eigen::MatrixXd A(3, 7);
        Eigen::VectorXd L(3);
        X.fill(0);
        X(0) = 1;
        Eigen::Vector3d x1_mass, x2_mass;
        x1_mass.fill(0);
        x2_mass.fill(0);
        for (auto pair : pairs) {
            x1_mass += pair.img_pt;
            x2_mass += pair.obj_pt;
        }
        x1_mass /= pairs.size();
        x2_mass /= pairs.size();

        int pt_count = 0;

        auto obj_fcn = [&](void*)->bool {
            Eigen::Vector3d Xhat, _Xhat;
            Xhat = pairs[pt_count].img_pt.transpose() - x1_mass;
            _Xhat = pairs[pt_count].obj_pt.transpose() - x2_mass;
            A.col(4) << 1, 0, 0;
            A.col(5) << 0, 1, 0;
            A.col(6) << 0, 0, 1;
            A.col(0) = Xhat.transpose();
            A.col(1) << -X(0) * Xhat(2), 0, X(0)* Xhat(0);
            A.col(2) << -Xhat(1) * sin(X(1)), Xhat(0)* sin(X(1)) - Xhat(2) * cos(X(1)), Xhat(1)* cos(X(1));
            A.col(2) *= X(0);
            A.col(3) << -Xhat(1) * cos(X(1)) * cos(X(2)) - Xhat(2) * sin(X(2)),
                Xhat(0)* cos(X(1))* cos(X(2)) + Xhat(2) * sin(X(1)) * cos(X(2)),
                Xhat(0)* sin(X(2)) - Xhat(1) * sin(X(1)) * cos(X(2));
            A.col(3) *= X(0);
            L = _Xhat - rotationMatrix(X(1), X(2), X(3)) * Xhat * X(0);

            if (++pt_count >= pairs.size()) {
                pt_count = 0;
                return true;
            }
            else {
                return false;
            }
        };

        opt::LeastSquareSolver solver(A, L, X);

        solver.setObjectFunction(obj_fcn);

        solver.method = method;

        solver.function_num = 3 * pairs.size();

        auto flag = solver.solve(nullptr);

        if (flag == opt::QuitFlag::MaxIterReached) {
            std::cout << "ERROR: Max iteration reached. Result unfeasable." << std::endl;
            return Eigen::VectorXd(3);
        }

        rotation = rotationMatrix(X(1), X(2), X(3));
        scale = X(0);
        translation = X.tail(3) + x2_mass - X(0) * rotationMatrix(PhyOmegaKappa, X(1), X(2), X(3)) * x1_mass;
    }
    Eigen::Vector3d error;
    error.fill(0);
    
    for (auto pair : pairs) {
        auto e = pair.obj_pt.transpose()-scale*rotation*pair.img_pt.transpose()-translation;
        error(0) += e(0) * e(0);
        error(1) += e(1) * e(1);
        error(2) += e(2) * e(2);
    }
    error /= pairs.size();
    error << sqrt(error(0)), sqrt(error(1)), sqrt(error(2));
    return error;
}

Eigen::VectorXd photo::sparseBundleAdjustment(cameraParas& camera, std::vector<std::map<int, Pair<2, 3>>>& pts, int num_img_pts, int num_obj_pts, std::vector<Eigen::VectorXd>& extrinsic_elements, opt::Method method)
{
    int num_photos = pts.size();
    int x0 = camera.cx();
    int y0 = camera.cy();
    double f = camera.fx();
    Eigen::MatrixXd _A(2, 6 * num_photos + 3 * num_obj_pts);
    Eigen::VectorXd _L(2);
    Eigen::VectorXd _X(6 * num_photos + 3 * num_obj_pts);

    for (int i = 0;i < num_photos;i++) {
        _X.block(i * 6, 0, 6, 1) = extrinsic_elements[i];
        for (auto pt : pts[i]) {
            _X.block(6 * num_photos + pt.first * 3, 0, 3, 1) = pt.second.obj_pt;
        }
    }

    int photo_id = 0, n = 0;
    auto pt = pts[photo_id].begin();
    auto obj_fcn = [&](void*)->bool {
        Eigen::Matrix3d R = rotationMatrix(_X(photo_id * 6 + 3), _X(photo_id * 6 + 4), _X(photo_id * 6 + 5));
        double Xslash = (pt->second.obj_pt - _X.block(photo_id * 6, 0, 3, 1).transpose()) * R.col(0);
        double Yslash = (pt->second.obj_pt - _X.block(photo_id * 6, 0, 3, 1).transpose()) * R.col(1);
        double Zslash = (pt->second.obj_pt - _X.block(photo_id * 6, 0, 3, 1).transpose()) * R.col(2);
        double x = pt->second.img_pt[0];
        double y = pt->second.img_pt[1];
        _L(0) = x - x0 + f * Xslash / Zslash;
        _L(1) = y - y0 + f * Yslash / Zslash;
        for (int i = 0;i < 3;i++) {
            _A(0, i + 6 * photo_id) = (R(i, 0) * f + R(i, 2) * (x - x0)) / Zslash;
            _A(1, i + 6 * photo_id) = (R(i, 1) * f + R(i, 2) * (y - y0)) / Zslash;
        }
        _A(0, 3 + 6 * photo_id) = (y - y0) * sin(_X(photo_id * 6 + 4)) - ((x - x0) / f * ((x - x0) * cos(_X(photo_id * 6 + 5)) - (y - y0) * sin(_X(photo_id * 6 + 5))) + f * cos(_X(photo_id * 6 + 5))) * cos(_X(photo_id * 6 + 4));
        _A(1, 3 + 6 * photo_id) = -(x - x0) * sin(_X(photo_id * 6 + 4)) - ((y - y0) / f * ((x - x0) * cos(_X(photo_id * 6 + 5)) - (y - y0) * sin(_X(photo_id * 6 + 5))) - f * sin(_X(photo_id * 6 + 5))) * cos(_X(photo_id * 6 + 4));
        _A(0, 4 + 6 * photo_id) = -f * sin(_X(photo_id * 6 + 5)) - (x - x0) / f * ((x - x0) * sin(_X(photo_id * 6 + 5)) + (y - y0) * cos(_X(photo_id * 6 + 5)));
        _A(1, 4 + 6 * photo_id) = -f * cos(_X(photo_id * 6 + 5)) - (y - y0) / f * ((x - x0) * sin(_X(photo_id * 6 + 5)) + (y - y0) * cos(_X(photo_id * 6 + 5)));
        _A(0, 5 + 6 * photo_id) = y - y0;
        _A(1, 5 + 6 * photo_id) = -(x - x0);
        _A.block(0, 6 * num_photos + pt->first * 3, 2, 3) = -_A.block(0, 0, 2, 3);
        pt++;
        if (pt == pts[photo_id].end()) {
            photo_id++;
            pt = pts[photo_id].begin();
        }

        if (photo_id >= num_photos) {
            photo_id = 0;
            return true;
        }
        else {
            return false;
        }
    };
    opt::LeastSquareSolver solver(_A, _L, _X);

    solver.function_num = num_img_pts * 2;

    solver.setObjectFunction(obj_fcn);

    solver.method = method;

    auto flag = solver.solve(nullptr);

    if (flag == opt::QuitFlag::MaxIterReached) {
        std::cout << "ERROR: Max iteration reached. Result unfeasable." << std::endl;
        return Eigen::VectorXd(3);
    }
}

