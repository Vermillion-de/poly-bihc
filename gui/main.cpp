#include <cstdio>
#include <iostream>
#include "imgui_viewer.hpp"
#include "surface.hpp"
#include "coordinates.hpp"

int main(int argc, char** argv){
    static int __frame__ = 0;
    easy3d::init();
    easy3d::ViewerImGui viewer;
    viewer.set_background_color({1.0f, 1.0f, 1.0f, 1.0f});
    viewer.picked_point_size = 30;
    int cage_curve_width = 2;
    
    std::string dir, method="combine";
    if(argc == 2){
        dir = std::string(argv[1]);
        method="combine";
    } if(argc == 3){
        dir = std::string(argv[1]);
        method = std::string(argv[2]);
    } else{
        std::cout << "[ ERROR ] ./main [data_dir] [initial method]" << std::endl;
    }

    MeshLoader mesh;
    BezierLoader cage, cage0, _cage;

    auto init_data = [&](std::string dir){
        mesh = MeshLoader(dir+"/mesh.obj");
        cage = BezierLoader(dir+"/cage1.obj");
        cage0 = BezierLoader(dir+"/cage1.obj");
        _cage = BezierLoader(dir+"/cage0.obj");
    };
    init_data(dir);

    float c1[] = {177./255, 177./255, 255./255}; // color 1: purple in papers
    float c2[] = {10./255, 10./255, 10./255}; // color 2: dark gray in papers

    easy3d::vec4 mesh_color{c1[0], c1[1], c1[2], 1.0f};
    easy3d::vec4 cage_color{c2[0], c2[1], c2[2], 1.0f};
    easy3d::vec4 point_color{0.5f, 0.2f, 0.8f, 0.5f};

    bool mesh_show_smooth_shading = true;

    int filt = 5;

    struct coefi_v_n {
        // coefi for low order
        MatrixXd g1n, g1_, g2n, g2_;
        MatrixXd G1n, G1_, G2n, G2_;
        MatrixXd Ig, _G1n, _G1_;
        MatrixXd coefi_v_green,  coefi_n_green;
        MatrixXd coefi_v_bih,    coefi_n_bih;
        MatrixXd coefi_v_comb,   coefi_n_comb;
        MatrixXd coefi_v_bih_reg,coefi_n_bih_reg;
        std::string coord_type = "Combine";
        MatrixXd *coefi_v, *coefi_n;
        float n_r=1, v_r=1;
        VectorXd mat_n_r, mat_v_r;
        float comb_v = 1, comb_n = 1, comb_uni = 1;
        float reg_1 = 1, reg_2 = 0;

        // coefi for higher orders 
        MatrixXd mat_lm_to_vn;

        void update_comb(){
            coefi_v_comb=G1n + comb_v * _G1n;
            coefi_n_comb=G1_ + comb_n * _G1_;
        };

        void update_bih_reg(){
            const int c_v = g1n.cols(), c_n = g1_.cols();
            const int c_l = g2n.cols(), c_m = g2_.cols();
            const MatrixXd mv_0 = MatrixXd::Zero(g1n.rows(), g1n.cols()); 
            const MatrixXd mn_0 = MatrixXd::Zero(g1_.rows(), g1_.cols()); 
            mat_lm_to_vn = optimize({g1n-Ig,  g2n}, {g1_,   g2_}, {c_l, c_m}, 
                                    {mv_0, Ig-g1n}, {mn_0, -g1_}, {c_v, c_n}, {reg_1, pow(10, reg_2)});
            auto [_c_v, _c_n] = vnlm_to_vn({G1n, G1_, G2n, G2_}, mat_lm_to_vn);
            coefi_v_bih_reg = _c_v;
            coefi_n_bih_reg = _c_n;
            auto [d_v, d_n] = lm_to_vn({G1n, G1_}, mat_lm_to_vn);
        }

        void update_Gs(){
            auto [__G1n, __G1_] = method1__G1n__G1_(g1n, g1_, g2n, g2_, Ig, 
                                                    G1n, G1_, G2n, G2_);
            _G1n = __G1n, _G1_ = __G1_;
        }

        void update_vn_r(){
            int n_v = coefi_v_green.cols();
            int n_n = coefi_n_green.cols();
            mat_v_r = VectorXd::Constant(n_v, v_r);
            mat_n_r = VectorXd::Constant(n_n, n_r);
        }

        void init(std::string dir){
            load(dir+"/g/Ig.csv", Ig);
            load(dir+"/g/g1n.csv", g1n); load(dir+"/g/g1_.csv", g1_);
            load(dir+"/g/g2n.csv", g2n); load(dir+"/g/g2_.csv", g2_);
            load(dir+"/G/G1n.csv", G1n); load(dir+"/G/G1_.csv", G1_);
            load(dir+"/G/G2n.csv", G2n); load(dir+"/G/G2_.csv", G2_);
            load(dir+"/coefi_v_green.csv", coefi_v_green);
            load(dir+"/coefi_n_green.csv", coefi_n_green);
            load(dir+"/coefi_v_bih.csv", coefi_v_bih);
            load(dir+"/coefi_n_bih.csv", coefi_n_bih);
            update_Gs();
            update_comb();
            update_vn_r();
            update_bih_reg();
            if(coord_type == "combine"){
                coefi_v = &coefi_v_comb;
                coefi_n = &coefi_n_comb;
            }else if(coord_type == "Green"){
                coefi_v = &coefi_v_green;
                coefi_n = &coefi_n_green;
            }else if(coord_type == "bih"){
                coefi_v = &coefi_v_bih;
                coefi_n = &coefi_n_bih;
            }else if(coord_type == "bih_reg"){
                coefi_v = &coefi_v_bih_reg;
                coefi_n = &coefi_n_bih_reg;
            }
        };

        coefi_v_n(std::string dir, std::string method="combine"){
            coord_type = method;
            init(dir);
        }

    } coefi_mesh(dir+"/mesh", method);

    bool show_cage_mesh = false;
    bool show_cage_boundary = true;
    bool show_cage_vert = true;
    bool show_cage_frame = false;
    bool show_ctrl_points = true;
    bool show_mesh = true;
    bool should_wait_others = false;
    int cage_resolution = 50;
    int fps = 30;
    real3 cur_point{0, 0, 0};
    MatrixXd cur_rot = MatrixXd::Identity(3, 3);
    std::string cur_add_pt_type = "pos";

    auto viewer_add_mesh = [&](std::string name, MeshLoader m, easy3d::vec4 color){
        if(show_mesh){
            viewer.add_mesh(name, 
                [&](int idx, int i){ 
                    // if(i == 0) return filter(m.vs(idx, i), 1.78); 
                    // if(i == 1) return filter(m.vs(idx, i), 4.20); 
                    // if(i == 2) return filter(m.vs(idx, i), 1.78); 
                    return filter(m.vs(idx, i), 5);
                }, m.vs.rows(), 
                [&](int idx, int i){ return m.fs[idx][i]; }, m.fs.size(),
                color, mesh_show_smooth_shading
            );
        }
        else {
            viewer.add_mesh(name, [&](int, int){return 0;}, 0, [](int, int){return 0;}, 0);
        }
        viewer.set_manager_bool(name+".edge_visible", false);
        viewer.set_manager_bool(name+".transparency", false);
    };

    auto viewer_add_cage_frame = [&](std::string name, BezierLoader m, easy3d::vec4 color){
        // std::cout << "Add Graph: " << name << std::endl; 
        MeshLoader mesh = m.to_skeleton();
        vector<vector<double>> none{};
        viewer.add_points(name+"_vert", show_cage_vert? to_double(m.vs) : none, 20, color);
        // viewer_add_mesh(name, mesh, color, false);
        if(show_cage_frame){
            auto [v_graph, e_graph] = m.to_graph();
            viewer.add_graph(name + "_frame", v_graph, e_graph, 3, 20, color);
        }
        else viewer.add_graph(name + "_frame", {}, {}, 3, 20, color);
        viewer.set_manager_bool(name+"_frame.edge_visible", true);
    };

    auto viewer_add_cage_mesh = [&](std::string name, BezierLoader m, easy3d::vec4 color){
        std::cout << "Add Cage: " << name << std::endl; 
        name += "_mesh";
        auto mesh = m.to_mesh(show_cage_mesh? cage_resolution : 0);
        viewer_add_mesh(name, mesh, {1.f, 1.f, 1.f, 1.f});
        viewer.set_manager_bool(name+".edge_visible", false);
        viewer.set_manager_bool(name+".transparency", true);
        auto [v_bound, e_bound] = m.to_skeleton(show_cage_boundary? cage_resolution : 0);
        viewer.add_graph(name + "_boundary", v_bound, e_bound, cage_curve_width, 3, color, cage_color);
        viewer.set_manager_bool(name+"_boundary.edge_visible", true);
    };


    std::string interface_type = "Vertex";

    auto vertex_update = [&](){
        mesh.vs = *(coefi_mesh.coefi_v) * rowwise_prod(cage.vs_mat, coefi_mesh.mat_v_r) + 
                  *(coefi_mesh.coefi_n) * rowwise_prod(cage.ns_mat, coefi_mesh.mat_n_r);
    };

    // bool first_use = false;
    // viewer.manage_checkbox("First_Use", &first_use, [&](){});
    auto update = [&](){
        // if(first_use) 
        if(interface_type == "None"){

        }else if(interface_type == "Vertex") 
            vertex_update();
        if(interface_type != "None" )
            viewer_add_mesh("Mesh", mesh, mesh_color);
        viewer_add_cage_mesh("Cage", cage, cage_color);
        viewer_add_cage_frame("Cage", cage, cage_color);
    };


    // Call Back:
    viewer.manage_picked_point_cloud("Cage_vert", [&](vector<int> ids, vector<vector<double>> pos, vector<vector<double>> rot){
        // cur_point = {pos[0], pos[1], pos[2]};
        if(ids.size() == 1) {
            cur_point = {pos[0][0], pos[0][1], pos[0][2]};
            std::cout << "Select ID: " << ids[0] << ", at (" << pos[0][0] << ", " << pos[0][1] << ", " << pos[0][2] << ")" << std::endl;
        };
        cage.update_vert(ids, pos);
        for(int r = 0; r < 3; r++) for(int c = 0; c < 3; c++)
            cur_rot(r, c) = rot[r][c];
        update();
    });

    viewer.manage_int("Cage Boundary Width", &cage_curve_width, 0, 10, [&](){
        viewer_add_cage_mesh("Cage", cage, cage_color);
    });

    string tgt_id;
    vector<string> saved_tgts = list_dir(dir+"/saved"); 
    // if(saved_tgts.size() != 0) tgt_id = saved_tgts[0];

    viewer.manage_button("Save_Screen", [&](){
        const std::string &dir = "../screenshot/data/" + base_name(string(argv[1]));
        auto id = utils::prefix();
        auto folder = dir + "/" + id; 
        std::filesystem::create_directories(folder);
        cage.save(folder+"/cage1_bzp.obj");
        cage0.save(folder+"/cage0_bzp.obj");
        cage.to_graph(folder+"/cage1.obj"); // save target cage frame
        cage0.to_graph(folder+"/cage0.obj"); // save target cage frame
        cage.to_skeleton(folder+"/cage1_bound.obj"); // save target cage lines 
        cage0.to_skeleton(folder+"/cage0_bound.obj"); // save target cage lines
        cage.to_mesh(cage_resolution).save(folder+"/cage1_mesh.obj"); // save target cage
        cage0.to_mesh(cage_resolution).save(folder+"/cage0_mesh.obj"); // save initial cage
        MeshLoader(dir+"/mesh.obj").save(folder+"/mesh0.obj");
        mesh.vs = coefi_mesh.coefi_v_green * rowwise_prod(cage.vs_mat, coefi_mesh.mat_v_r) + // result of Green 
                  coefi_mesh.coefi_n_green * rowwise_prod(cage.ns_mat, coefi_mesh.mat_n_r);
        mesh.save(folder+"/mesh_green.obj");
        mesh.vs = coefi_mesh.coefi_v_bih * rowwise_prod(cage.vs_mat, coefi_mesh.mat_v_r) +  // result of Weber's BiH
                  coefi_mesh.coefi_n_bih * rowwise_prod(cage.ns_mat, coefi_mesh.mat_n_r);
        mesh.save(folder+"/mesh_bih.obj");
        mesh.vs = coefi_mesh.coefi_v_bih_reg * rowwise_prod(cage.vs_mat, coefi_mesh.mat_v_r) + // result of BiH's reg
                  coefi_mesh.coefi_n_bih_reg* rowwise_prod(cage.ns_mat, coefi_mesh.mat_n_r);
        mesh.save(folder+"/mesh_bih_reg.obj");
    });

    viewer.manage_button("Save", [&](){
        auto id = utils::prefix();
        auto folder = dir + "/saved/" + id; 
        std::filesystem::create_directories(folder);
        cage.save(folder+"/cage_tmp.obj");
        cage.to_mesh(cage_resolution).save(folder+"/cage_mesh_tmp.obj");
        viewer.save_camera_parameter(folder+"/camera.txt");
        mesh.save(folder+"/mesh_tmp.obj");
        saved_tgts.push_back(id);
        viewer.manage_selectbox("Saved Targeds", saved_tgts, &tgt_id, [&](){
            auto folder = dir+"/saved/"+tgt_id;
            cage = BezierLoader(folder+"/cage_tmp.obj");
            update();
        });
    });

    viewer.manage_selectbox("Saved Targeds", saved_tgts, &tgt_id, [&](){
        auto folder = dir+"/saved/"+tgt_id;
        cage = BezierLoader(folder+"/cage_tmp.obj");
        update();
    });

    viewer.manage_button("Reset All", [&](){
        cage = cage0; 
        update();
    });


    viewer.manage_checkbox("Show Cage Mesh", &show_cage_mesh, [&](){
        viewer_add_cage_mesh("Cage", cage, cage_color); });

    viewer.manage_checkbox("Show Cage Boundary", &show_cage_boundary, [&](){
        viewer_add_cage_mesh("Cage", cage, cage_color); });

    viewer.manage_checkbox("Show Cage Vert", &show_cage_vert, [&](){
        viewer_add_cage_frame("Cage", cage, cage_color); });

    viewer.manage_checkbox("Show Cage Frame", &show_cage_frame, [&](){
        viewer_add_cage_frame("Cage", cage, cage_color); });

    viewer.manage_checkbox("Show Mesh", &show_mesh, [&](){ 
        viewer_add_mesh("Mesh", mesh, mesh_color);});

    viewer.manage_checkbox("Show Mesh Smooth Shading", &mesh_show_smooth_shading, [&](){
        viewer_add_mesh("Mesh", mesh, mesh_color);});

    int cage_degree = cage.get_degree();
    viewer.manage_int("Cage - Degree", &cage_degree, 1, 4, [&](){
        cage = _cage; 
        cage.change_degree(cage_degree);
        viewer_add_cage_mesh("Cage", cage, cage_color);
        viewer_add_cage_frame("Cage", cage, cage_color);
    });

    viewer.manage_int("Bezier Resolution", &cage_resolution, 0, 100, [&](){
        viewer_add_cage_mesh("Cage", cage, cage_color); });

    viewer.manage_float("Picked Point Size", &viewer.picked_point_size, 0, 200, [&](){ update();});

    std::string normal_type="uvw";
    viewer.manage_selectradio("Normal Type", {"uvw", "vert"}, &normal_type, [&](){
        if(normal_type == "uvw") cage.change_normal_type("uvw");
        if(normal_type == "vert") cage.change_normal_type("vert");
        update();
    });

    viewer.manage_selectradio("Interface Type", {"None", "Vertex"}, &interface_type, [&](){
        update();
    });

    viewer.manage_int("Optimize Filter (exp)", &filt, 0, 100, [&](){
        update();
    });

    auto manage_coefi = [&](coefi_v_n& c){
        viewer.manage_selectbox("Coordinate", {"Green", "BiH", "Combine", "BiH Reg"}, &c.coord_type, [&](){
            if(c.coord_type == "Green") 
                c.coefi_v = &c.coefi_v_green, c.coefi_n = &c.coefi_n_green;
            if(c.coord_type == "BiH")
                c.coefi_v = &c.coefi_v_bih, c.coefi_n = &c.coefi_n_bih;
            if(c.coord_type == "BiH Reg")
                c.coefi_v = &c.coefi_v_bih_reg, c.coefi_n = &c.coefi_n_bih_reg;
            if(c.coord_type == "Combine")
                c.coefi_v = &c.coefi_v_comb, c.coefi_n = &c.coefi_n_comb;
            update();
        });

        viewer.manage_float("v ratio", &c.v_r, -1, 3, [&](){ c.update_vn_r(); update(); });
        viewer.manage_float("n ratio", &c.n_r, -1, 3, [&](){ c.update_vn_r(); update(); });

        viewer.manage_float("v comb", &c.comb_v, -1, 3, [&](){ 
            c.coefi_v_comb =  c.G1n + c.comb_v * c._G1n; update(); });

        viewer.manage_float("n comb", &c.comb_n, -1, 3, [&](){
            c.coefi_n_comb =  c.G1_ + c.comb_n * c._G1_; update(); });

        viewer.manage_float("uni comb", &c.comb_uni, -1, 2, [&](){
        c.comb_v = c.comb_n = c.comb_uni; c.update_comb(); update(); });

        viewer.manage_float("reg lamb", &c.reg_2, -20, 20, [&](){
                        c.update_bih_reg();update(); });

        static string lamb = "0";
        std::vector<string> lambs;
        for(int i = -20; i <= 20; i += 2) lambs.push_back(std::to_string(i));
        viewer.manage_selectbox("Select reg lam", lambs, &lamb, [&](){
            c.reg_2 = std::stof(lamb);
            c.update_bih_reg(); update();
        });
    };
    manage_coefi(coefi_mesh);
    update();
    return viewer.run();
}
