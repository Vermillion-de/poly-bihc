#include "imgui_viewer.hpp"

#include <easy3d/core/model.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/graph.h>

#include <easy3d/renderer/camera.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/manipulator.h>
#include <easy3d/renderer/dual_depth_peeling.h>
#include <easy3d/renderer/manipulated_frame.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/text_renderer.h>
#include <easy3d/renderer/shape.h>

#include <easy3d/gui/picker_model.h>
#include <easy3d/gui/picker_surface_mesh.h>
#include <easy3d/gui/picker_point_cloud.h>
#include <easy3d/gui/picker.h>
#include <easy3d/util/file_system.h>
#include <easy3d/util/setting.h>

#define IMGUI_DEFINE_MATH_OPERATORS
#include <3rd_party/glfw/include/GLFW/glfw3.h> 
#include <3rd_party/imgui/misc/fonts/imgui_fonts_droid_sans.h>
#include <3rd_party/imgui/imgui.h>
#include <3rd_party/imgui/imgui_internal.h>
#include <3rd_party/imgui/backends/imgui_impl_glfw.h>
#include <3rd_party/imgui/backends/imgui_impl_opengl3.h>
#include <3rd_party/imgui/imgui_demo.cpp>

#include "implot.h"
#include "implot_internal.h"
#include "implot.cpp"
#include "implot_items.cpp"
#include "implot_demo.cpp"
// #include "imguizmo.h"

#include "manager.hpp"

#include <vector>
#include <queue>
#include <map>
#include <omp.h>

#define DEBUG_COUT {std::cout << "Fine Here" << __LINE__;}
using namespace std;

namespace easy3d{

ImGuiContext* ViewerImGui::context_ = nullptr;
ViewerImGui::ViewerImGui( const std::string& title, int samples ,
        int gl_major, int gl_minor, bool full_screen , bool resizable,
        int depth_bits, int stencil_bits, int width, int height
) : Viewer(title, samples, gl_major, gl_minor, full_screen, resizable, depth_bits, stencil_bits, width, height)
{
#if defined(_WIN32) && defined(_MSC_VER)
    glfwInit(); // Liangliang: the internal glfw won't be shared across dll boundaries (But seems ok on macOS. That is weird!)
#endif
}

ViewerImGui::~ViewerImGui() {
    ImGui_ImplOpenGL3_Shutdown(); ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext(context_);
    ImPlot::DestroyContext();
#if defined(_WIN32) && defined(_MSC_VER)
    glfwTerminate(); // Liangliang: the internal glfw won't be shared across dll boundaries (But seems ok on macOS. That is weird!)
#endif
}

void ViewerImGui::save_camera_parameter(std::string file){
    std::ofstream os(file);
    vec3 pos = camera()->position();
    vec3 up = camera()->upVector();
    vec3 dir = camera()->viewDirection();
    os << "Position = (" << pos.x << ", " << pos.y << ", " << pos.z << ")" << std::endl;
    os << "Up = (" << up.x << ", " << up.y << ", " << up.z << ")" << std::endl;
    os << "Dir = (" << dir.x << ", " << dir.y << ", " << dir.z << ")" << std::endl;
    os.close();
}

static DualDepthPeeling* __transparancy_render__;
void ViewerImGui::init() {
    Viewer::init();
    if (!context_) {
        IMGUI_CHECKVERSION();
        context_ = ImGui::CreateContext();
        ImPlot::CreateContext();
        const char* glsl_version = "#version 150";
        ImGui_ImplGlfw_InitForOpenGL(window_, true);
        ImGui_ImplOpenGL3_Init(glsl_version);
        ImGuiIO& io = ImGui::GetIO();
        io.WantCaptureKeyboard = true;
        io.WantTextInput = true;
        io.IniFilename = nullptr;
        ImGui::StyleColorsLight();
        ImGuiStyle& style = ImGui::GetStyle();
        style.FrameRounding = 5.0f;
        style.Alpha = 0.9f;
        reload_font();
        __transparancy_render__ = new DualDepthPeeling(camera());
    }
    if (!drawable_axes_) {
        float base = 0.25f;   // the cylinder length, relative to the allowed region
        float head = 0.1f;   // the cone length, relative to the allowed region
        std::vector<vec3> points, normals, colors;
        shape::create_cylinder(0.01, 10, vec3(0, 0, 0), vec3(base, 0, 0), vec3(1, 0, 0), points, normals, colors);
        shape::create_cylinder(0.01, 10, vec3(0, 0, 0), vec3(0, base, 0), vec3(0, 1, 0), points, normals, colors);
        shape::create_cylinder(0.01, 10, vec3(0, 0, 0), vec3(0, 0, base), vec3(0, 0, 1), points, normals, colors);
        shape::create_cone(0.02, 20, vec3(base, 0, 0), vec3(base + head, 0, 0), vec3(1, 0, 0), points, normals,
                                colors);
        shape::create_cone(0.02, 20, vec3(0, base, 0), vec3(0, base + head, 0), vec3(0, 1, 0), points, normals,
                                colors);
        shape::create_cone(0.02, 20, vec3(0, 0, base), vec3(0, 0, base + head), vec3(0, 0, 1), points, normals,
                                colors);
        shape::create_sphere(vec3(0, 0, 0), 0.02, 20, 20, vec3(0, 1, 1), points, normals, colors);
        const_cast<ViewerImGui*>(this)->drawable_axes_ = new TrianglesDrawable("corner_axes");
        drawable_axes_->update_vertex_buffer(points);
        drawable_axes_->update_normal_buffer(normals);
        drawable_axes_->update_color_buffer(colors);
        drawable_axes_->set_property_coloring(State::VERTEX);
    }
}

void ViewerImGui::reload_font(int font_size)
{
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();
    io.Fonts->AddFontFromMemoryCompressedTTF(droid_sans_compressed_data, droid_sans_compressed_size, static_cast<float>(font_size) * dpi_scaling());
    io.FontGlobalScale = width() * 1.0f / framebuffer_width();
    ImGui_ImplOpenGL3_DestroyDeviceObjects();
}

void ViewerImGui::post_resize(int w, int h) {
    Viewer::post_resize(w, h);
    if (context_) {
        ImGui::GetIO().DisplaySize.x = float(w);
        ImGui::GetIO().DisplaySize.y = float(h);
    }
}

void ViewerImGui::pre_draw() {
    // Calculate ms/frame
    static double last_time = 0;
    static int frame_counter = 0;
    double current_time = glfwGetTime();
    ++frame_counter;
    if (current_time - last_time >= 1.0f) {
        refresh_rate = double(frame_counter) / (current_time - last_time);
        frame_counter = 0;
        last_time = current_time;
    }

    ImGui_ImplGlfw_NewFrame();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui::NewFrame();   
    Viewer::pre_draw(); 
}

void init(){ 
    initialize(); 
}

// callbacks.
bool ViewerImGui::callback_event_cursor_pos(double x, double y) {
    if (ImGui::GetIO().WantCaptureMouse) return true;
    else return Viewer::callback_event_cursor_pos(x, y);
}

bool ViewerImGui::callback_event_mouse_button(int button, int action, int modifiers) {
    // std::cout << "Mouse Pressed CallBack!" << std::endl;
    if (ImGui::GetIO().WantCaptureMouse) return true;
    else return Viewer::callback_event_mouse_button(button, action, modifiers);
}

bool ViewerImGui::callback_event_keyboard(int key, int action, int modifiers) {
    // if(custom_keycallback(key, action, modifiers)) return true;
    if (ImGui::GetIO().WantCaptureKeyboard) return true;
    else return Viewer::callback_event_keyboard(key, action, modifiers);
}

bool ViewerImGui::callback_event_character(unsigned int codepoint) {
    if (ImGui::GetIO().WantCaptureKeyboard) return true;
    else return Viewer::callback_event_character(codepoint);
}

bool ViewerImGui::callback_event_scroll(double dx, double dy) {
    if (ImGui::GetIO().WantCaptureMouse) return true;
    else return Viewer::callback_event_scroll(dx, dy);
}

///////////////////////// More than customized /////////////////////



static std::vector<std::string> _mng_log;
static vector<string> _mng_sec; // ImGui Sections
static sec_manager<function<void()>> _mng_para_callback;
static std::map<std::string, ImVec4> manager_color;
static std::map<std::string, bool> manager_bool;
static std::map<std::string, float> manager_system_float;
static map<string, PointCloud*> manager_PointCloud;
static map<string, Graph*> manager_Graph;
static map<string, SurfaceMesh*> manager_SurfaceMesh;


static sec_manager<tuple<float*, float, float>> _mng_para_float;
void ViewerImGui::manage_float(string words, float* value, float f_min, float f_max, function<void()> callback){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_para_float.add_object(sec, name, {value, f_min, f_max});
    _mng_para_callback.add_object(sec, name, callback);
}

static sec_manager<tuple<int*, int, int>> _mng_para_int;
void ViewerImGui::manage_int(string words, int* value, int i_min, int i_max, function<void()> callback){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_para_int.add_object(sec, name, {value, i_min, i_max});
    _mng_para_callback.add_object(sec, name, callback);
}

static sec_manager<tuple<int*, int, int>> _mng_para_int2;
void ViewerImGui::manage_int2(string words, int* value, int i_min, int i_max, function<void()> callback){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_para_int2.add_object(sec, name, {value, i_min, i_max});
    _mng_para_callback.add_object(sec, name, callback);
}

static sec_manager<function<void()>> _mng_button;
void ViewerImGui::manage_button(string words, function<void()> action){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_button.add_object(sec, name, action);
}

static sec_manager<tuple<bool*,function<void()>>> _mng_checkbox;
void ViewerImGui::manage_checkbox(string words, bool* val, function<void()> action){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_checkbox.add_object(sec, name, {val, action});
}

static map<string, function<void(int, vector<double>, vector<vector<double>>)>> _mng_mani_0;
void ViewerImGui::manage_picked_faced_point(string name, function<void(int, vector<double>, vector<vector<double>>)> action){
    _mng_mani_0[name] = action;
    set_manager_bool(name+".pick", true);
}

static map<string, function<void(vector<int>, vector<vector<double>>, vector<vector<double>>)>> _mng_mani_1;
void ViewerImGui::manage_picked_point_cloud(string name, function<void(vector<int>, vector<vector<double>>, vector<vector<double>>)> action){
    _mng_mani_1[name] = action;
    set_manager_bool(name+".pick", true);
}

static sec_manager<tuple<string*, vector<string>>> _mng_selectbox;
void ViewerImGui::manage_selectbox(string words, vector<string> items, string* val, function<void()> callback){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_selectbox.add_object(sec, name, {val, items});
    _mng_para_callback.add_object(sec, name, {callback});
}

static sec_manager<tuple<string*, vector<string>>> _mng_selectradio;
void ViewerImGui::manage_selectradio(string words, vector<string> items, string* val, function<void()> callback){
    auto parts = split(words, '.');
    auto sec = parts[0], name = parts[1];
    unique_append(_mng_sec, sec);
    _mng_selectradio.add_object(sec, name, {val, items});
    _mng_para_callback.add_object(sec, name, {callback});
}

// on keycallback
static map<char, function<void()>> key_callback;
void ViewerImGui::manage_keycallback(char val, function<void()> action){
    key_callback[val] = action;   
}

void ViewerImGui::manage_clear(std::string type){
    if (type == "int") _mng_para_int.clear(); 
    else if(type == "float") _mng_para_float.clear(); 
    else if(type == "selectbox") _mng_selectbox.clear();
    else if(type == "button") _mng_button.clear();
    else if(type == "callback") _mng_para_callback.clear();
}

void ViewerImGui::manage_clear(std::string type, std::string sec, std::string name){
    if (type == "int")  _mng_para_int.clear(sec, name); 
    else if(type == "float") _mng_para_float.clear(sec, name); 
    else if(type == "selectbox") _mng_selectbox.clear(sec, name);
    else if(type == "button") _mng_button.clear(sec, name);
    else if(type == "callback") _mng_para_callback.clear(sec, name);
}

void ViewerImGui::custom_imgui(){
    static bool open_main_window_ptr = true;
    auto cond = ImGuiCond_FirstUseEver; // ImGuiCond_Always
    ImGui::SetNextWindowPos(ImVec2(width() * 0, height() * 0.03333), cond, ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(width() * 0.300, height() * 0.950), cond);
    ImGui::Begin("Custom ImGui", &open_main_window_ptr);
    for(auto sec : _mng_sec){
        ImGui::SeparatorText(sec.c_str());

        _mng_checkbox.apply(sec, [&](string name, tuple<bool*, function<void()>> data){
            auto val = get<0>(data); auto action = get<1>(data);
            if(ImGui::Checkbox(name.c_str(), val)) 
                action();
        });

        _mng_button.apply(sec, [&](string name, function<void()> action){
            if(ImGui::Button(name.c_str())) 
                action();
        });

        _mng_selectradio.apply(sec, [&](string name, tuple<string*, vector<string>> select_box){
            auto [tgt, items] = select_box;
            ImGui::Text(name.c_str());
            for(auto it : items){
                ImGui::SameLine();
                if(ImGui::RadioButton(it.c_str(), it == *tgt)){
                    *tgt = it;
                    _mng_para_callback.apply(sec, name, [](function<void()> f){f();});
                }
            }
        });

        _mng_selectbox.apply(sec, [&](string name, tuple<string*, vector<string>> select_box){
            auto tgt =   std::get<0>(select_box);
            auto items = std::get<1>(select_box);
            if(ImGui::BeginCombo(name.c_str(), (*tgt).c_str())){
                for(auto sel_item : items){
                    bool is_selected = (sel_item == *tgt);
                    if(ImGui::Selectable(sel_item.c_str(), is_selected))
                    {
                        *tgt = sel_item;
                        _mng_para_callback.apply(sec, name, [](function<void()> f){f();});
                    }
                    if(is_selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        });

        _mng_para_int2.apply(sec, [&](string name, tuple<int*, int, int> data){
            auto val = std::get<0>(data); auto imin = std::get<1>(data); auto imax = std::get<2>(data);
            if(ImGui::SliderInt2(name.c_str(), val, imin, imax))
                _mng_para_callback.apply(sec, name, [](function<void()> f){f();});
        });

        _mng_para_float.apply(sec, [&](string name, tuple<float*, float, float> data){
            auto val = std::get<0>(data); auto fmin = std::get<1>(data); auto fmax = std::get<2>(data);
            if(ImGui::SliderFloat(name.c_str(), val, fmin, fmax))
                _mng_para_callback.apply(sec, name, [](function<void()> f){f();});
        });

        _mng_para_int.apply(sec, [&](string name, tuple<int*, int, int> data){
            auto val = std::get<0>(data); auto imin = std::get<1>(data); auto imax = std::get<2>(data);
            if(ImGui::SliderInt(name.c_str(), val, imin, imax))
                _mng_para_callback.apply(sec, name, [](function<void()> f){f();});
        });
    }

    ImGui::SeparatorText("Logs");
    std::string message; for(auto log : _mng_log) message += log + "\n";
    ImGui::InputTextMultiline("##logs", &message[0], message.size(), ImVec2(-FLT_MIN, ImGui::GetTextLineHeight() * 7), ImGuiInputTextFlags_ReadOnly);
    refresh_rate = int(refresh_rate * 100) * 1. / 100;
    auto fps = "FPS:" + std::to_string(refresh_rate);
    ImGui::SeparatorText(fps.c_str());
    render_settings();
    ImGui::End();
}

void ViewerImGui::render_settings(){
    static bool alpha_preview = true, alpha_half_preview = true, drag_and_drop = true;
    static bool show_system_imgui = false;

    ImGuiColorEditFlags misc_flags = (drag_and_drop ? 0 : ImGuiColorEditFlags_NoDragDrop) |
                                     (alpha_half_preview ? ImGuiColorEditFlags_AlphaPreviewHalf : (alpha_preview ? ImGuiColorEditFlags_AlphaPreview : 0)) |
                                     ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel;

    ImGui::SetNextItemOpen(true, ImGuiCond_FirstUseEver);
    if(ImGui::TreeNodeEx("Rendering Settings")){
        ImGui::SeparatorText("Surfaces");
        for(auto [name, mesh] : manager_SurfaceMesh){
            auto edges = mesh->renderer()->get_lines_drawable("edges");
            auto faces = mesh->renderer()->get_triangles_drawable("faces");

            ImGui::SetNextItemOpen(true, ImGuiCond_FirstUseEver);
            if(ImGui::TreeNodeEx(name.c_str())){
                auto is_equal = [](ImVec4 a, vec4 b) -> bool {
                    if(a.x != b.x || a.y != b.y) return false;
                    if(a.z != b.z || a.w != b.w) return false;
                    return true;
                };

                auto name_edge_transparency = name + ".transparency";
                if(ImGui::Checkbox("transparency", &(manager_bool[name_edge_transparency]))){
                }
                ImGui::SameLine();

                manager_bool[name] = mesh->renderer()->is_visible();
                if(ImGui::Checkbox("Visible", &(manager_bool[name])))
                    mesh->renderer()->set_visible(manager_bool[name]);
                ImGui::SameLine();

                auto c = faces->color();
                manager_color[name] = {c.x, c.y, c.z, c.w};
                if(ImGui::ColorEdit4("Color", (float*)&manager_color[name], misc_flags)){
                    ImVec4 c = manager_color[name];
                    faces->set_uniform_coloring({c.x, c.y, c.z, c.w});
                }
                ImGui::SameLine();

                ImGui::Text("Color");
                ImGui::SameLine();

                auto name_edge_visble = name + ".edge_visible";
                manager_bool[name_edge_visble] = edges->is_visible();
                if(ImGui::Checkbox("Edges", &(manager_bool[name_edge_visble])))
                    edges->set_visible(manager_bool[name_edge_visble]);
                ImGui::SameLine();

                auto name_back_color = name + ".back_color";
                manager_bool[name_back_color] = faces->distinct_back_color();
                if(ImGui::Checkbox("Back Color", &(manager_bool[name_back_color]))){
                    faces->set_distinct_back_color(manager_bool[name_back_color]);
                    auto c = faces->back_color();
                }

                if(manager_bool[name_back_color]){
                    auto c = faces->back_color();
                    manager_color[name_back_color] = {c.x, c.y, c.z, c.w}; 
                    if(ImGui::ColorEdit4("Back Color", (float*)&manager_color[name_back_color], misc_flags)){
                        ImVec4 c = manager_color[name_back_color];
                        faces->set_back_color({c.x, c.y, c.z, c.w});
                    }
                    ImGui::SameLine();
                    ImGui::Text("Color");
                }
                ImGui::TreePop();
            }
        }

        ImGui::SeparatorText("Graphs");
        for(auto [name, graph] : manager_Graph){
            auto ptr_l = graph->renderer()->get_lines_drawable("edges");
            auto ptr_p = graph->renderer()->get_points_drawable("vertices");

            ImGui::SetNextItemOpen(true, ImGuiCond_FirstUseEver);
            if(ImGui::TreeNodeEx(name.c_str())){
                
                ImGui::SeparatorText("Edges");

                manager_system_float[name] = ptr_l->line_width();
                if(ImGui::SliderFloat("Width", &manager_system_float[name], 0, 50)){
                    float width = manager_system_float[name];
                    ptr_l->set_line_width(width);
                }

                manager_bool[name] = ptr_l->is_visible();
                if(ImGui::Checkbox("Visible", &(manager_bool[name])))
                    ptr_l->set_visible(manager_bool[name]);
                ImGui::SameLine();

                auto c = ptr_l->color();
                manager_color[name] = {c.r, c.g, c.b, c.a};
                if(ImGui::ColorEdit4("Color", (float*)&manager_color[name], misc_flags)){
                    ImVec4 c = manager_color[name];
                    ptr_l->set_uniform_coloring({c.x, c.y, c.z, c.w});
                }
                ImGui::SameLine();
                ImGui::Text("Color");

                ImGui::SeparatorText("Points");

                manager_system_float[name] = ptr_p->point_size();
                if(ImGui::SliderFloat("Size", &manager_system_float[name], 0, 50)){
                    float size = manager_system_float[name];
                    ptr_p->set_point_size(size);
                }

                manager_bool[name] = ptr_p->is_visible();
                if(ImGui::Checkbox("Visible", &(manager_bool[name])))
                    ptr_p->set_visible(manager_bool[name]);
                ImGui::SameLine();

                auto c_ = ptr_p->color();
                manager_color[name] = {c_.r, c_.g, c_.b, c_.a};
                if(ImGui::ColorEdit4("Color", (float*)&manager_color[name], misc_flags)){
                    ImVec4 c = manager_color[name];
                    ptr_p->set_uniform_coloring({c.x, c.y, c.z, c.w});
                }
                ImGui::SameLine();
                ImGui::Text("Color");

                ImGui::TreePop();
            }

        }

        ImGui::SeparatorText("Color Options");
        if(ImGui::TreeNodeEx("Color Options")){
            ImGui::Checkbox("With Alpha Preview", &alpha_preview);
            ImGui::Checkbox("With Half Alpha Preview", &alpha_half_preview);
            ImGui::Checkbox("With Drag and Drop", &drag_and_drop);
            ImGui::TreePop();
        }

        auto c = this->background_color();
        manager_color["Background Color"] = {c.r, c.g, c.b, c.a};
        if(ImGui::ColorEdit4("Background Color", (float*)&manager_color["Background Color"], misc_flags)){
            ImVec4 c = manager_color["Background Color"];
            this->set_background_color({c.x, c.y, c.z, c.w});
        }
        ImGui::SameLine();
        ImGui::Text("Background Color");
        ImGui::SameLine();

        if(ImGui::Button("Centering Object")){
            this->fit_screen();
            // camera()->setPosition(from)
            // camera()->centerScene();
        }

        ImGui::TreePop();
    }
}

static map<string, pair<vector<vector<double>>, vector<vector<double>>> > _mng_plot_2d;
void ViewerImGui::add_implot_lines(std::string name, std::vector<std::vector<std::vector<double>>> points, float line_width) 
{
    vector<vector<double>> plot_x(points.size());
    vector<vector<double>> plot_y(points.size());
    for(int i = 0; i < points.size(); i++){
        plot_x[i].resize(points[i].size());
        plot_y[i].resize(points[i].size());
        for(int j = 0; j < points[i].size(); j++){
            plot_x[i][j] = points[i][j][0];
            plot_y[i][j] = points[i][j][1];
        }
    }
    _mng_plot_2d[name] = {plot_x, plot_y};
}

void ViewerImGui::log_imgui(std::string log){
    if(_mng_log.size() == 5){
        std::vector<std::string> mng_log;
        for(int i = 1; i < 5; i++)
            mng_log.push_back(_mng_log[i]);
        mng_log.push_back(log);
        _mng_log = mng_log;
    }
    else {
        _mng_log.push_back(log);
    }
}

void ViewerImGui::add_handle_surf( std::string name, std::function<vector<double>(double, double)> pos,
                                    vector<double> u01, vector<double> v01, int res, vec4 color){
    double u0 = u01[0], u1 = u01[1];
    double v0 = v01[0], v1 = v01[1];
    
    SurfaceMesh* mesh; 
    
    if(manager_SurfaceMesh.find(name) != manager_SurfaceMesh.end()){
        mesh = manager_SurfaceMesh[name];
        delete_model(mesh);
    }

    mesh = new SurfaceMesh; mesh->set_name(name);
    manager_SurfaceMesh[name] = mesh;

    vector<SurfaceMesh::Vertex> vert;
    for(int i=0; i < res + 1; i++) for(int j=0; j < res + 1; j++) {
        auto p = pos(u0 + (u1 - u0) / res * i, v0 + (v1 - v0) / res * j);
        auto v = mesh->add_vertex({float(p[0]), float(p[1]),float(p[2])});
        vert.push_back(v);
    }
    auto idx = [&](int i, int j){return vert[i * (res + 1) + j];};
    for(int i=0; i < res; i++) for(int j=0; j < res; j++){ 
        mesh->add_triangle(idx(i, j), idx(i, j+1), idx(i+1, j));
        mesh->add_triangle(idx(i, j+1), idx(i+1, j+1), idx(i+1, j));
    } 
    auto model = add_model(mesh, true);
    auto faces = model->renderer()->get_triangles_drawable("faces");
    faces->set_uniform_coloring(color);
    faces->set_distinct_back_color(false);
}

void ViewerImGui::add_mesh(std::string name, std::function<float(int, int)> fv, int nv,
                                std::function<int(int, int)> ff, int nf, vec4 color, bool smooth_shading){
    SurfaceMesh* mesh; 
    if(manager_SurfaceMesh.find(name) != manager_SurfaceMesh.end()){
        mesh = manager_SurfaceMesh[name];
        delete_model(mesh);
    }

    mesh = new SurfaceMesh; mesh->set_name(name);
    manager_SurfaceMesh[name] = mesh;

    vector<SurfaceMesh::Vertex> vert;
    for(int i = 0; i < nv; i++){
        auto v = mesh->add_vertex({fv(i, 0), fv(i, 1), fv(i, 2)});
        vert.push_back(v);
    }

    for(int i = 0; i < nf; i++){
        mesh->add_triangle(vert[ff(i, 0)], vert[ff(i, 1)], vert[ff(i, 2)]);
    }
    auto model = add_model(mesh, true);
    auto faces = model->renderer()->get_triangles_drawable("faces");
    faces->set_uniform_coloring(color);
    faces->set_distinct_back_color(false);
    faces->set_smooth_shading(smooth_shading);
}

void ViewerImGui::add_graph(std::string name, std::vector<std::vector<double>> pts, std::vector<std::vector<int>> edges, float line_width, float point_size, vec4 color_point, vec4 color_edge){
    Graph* graph;
    if(manager_Graph.find(name) != manager_Graph.end()){
        graph = manager_Graph[name];
        delete_model(graph);
    }

    graph = new Graph; graph->set_name(name);
    manager_Graph[name] = graph;

    std::vector<easy3d::Graph::Vertex> vert(pts.size());
    for(int i = 0; i < pts.size(); i++)
        vert[i] = graph->add_vertex(vec3{float(pts[i][0]), float(pts[i][1]), float(pts[i][2])});
    for(auto edge : edges) for(int i = 0; i < edge.size()-1; i++)
        graph->add_edge(vert[edge[i]], vert[edge[i+1]]);

    add_model(graph);
    auto ld = graph->renderer()->get_lines_drawable("edges");
    ld->set_line_width(line_width);
    ld->set_impostor_type(LinesDrawable::CYLINDER);
    ld->set_uniform_coloring(color_edge);
    auto pd = graph->renderer()->get_points_drawable("vertices");
    pd->set_point_size(point_size);
    pd->set_uniform_coloring(color_point);
    pd->set_impostor_type(PointsDrawable::SPHERE);
}

void ViewerImGui::add_points(std::string name, std::vector<std::vector<double>> pts, float point_size, vec4 color){
    PointCloud* ptc;
    if(manager_PointCloud.find(name) != manager_PointCloud.end()){
        ptc = manager_PointCloud[name];
        delete_model(ptc);
    }
    ptc = new PointCloud; ptc->set_name(name);
    manager_PointCloud[name] = ptc;
    for(auto pt : pts) 
        ptc->add_vertex({float(pt[0]), float(pt[1]), float(pt[2])});
    add_model(ptc);
    auto pd = ptc->renderer()->get_points_drawable("vertices");
    pd->set_uniform_coloring(color); pd->set_point_size(point_size);
    pd->set_impostor_type(PointsDrawable::SPHERE);
}

map<string, mat4> manager_axes_wordMat;
void ViewerImGui::set_axes(std::string name, vector<vector<double>> wordMat){
    // std::cout << "Fine Here " << __LINE__ << std::endl;
    if(wordMat.size() != 0){
        auto m = wordMat;
        // std::cout << "Fine Here " << __LINE__ << std::endl;
        manager_axes_wordMat[name] = {
            float(m[0][0]), float(m[0][1]), float(m[0][2]), float(m[0][3]), 
            float(m[1][0]), float(m[1][1]), float(m[1][2]), float(m[1][3]), 
            float(m[2][0]), float(m[2][1]), float(m[2][2]), float(m[2][3]), 
            float(m[3][0]), float(m[3][1]), float(m[3][2]), float(m[3][3]), 
        };
        // std::cout << "Fine Here " << __LINE__ << std::endl;
    }
    else {
        if(manager_axes_wordMat.find(name) != manager_axes_wordMat.end())
            manager_axes_wordMat.erase(name);
    }
}

easy3d::Model* _picked_pos=nullptr;
std::string picking_mode = "faced_point";

// picking for mesh 
easy3d::SurfaceMesh* _picked_mesh=nullptr;
easy3d::SurfaceMesh::Vertex _picked_vert;

// picking for point cloud
std::string _picked_point_cloud_name;
vector<int> _picked_points_idx;
vector<vector<double>> _picked_points_pos;
easy3d::Polygon2 _picked_polygon;

vector<vector<double>> _I_mat4_ = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
}; 

bool ViewerImGui::mouse_press_event(int x, int y, int button, int modifiers) {
    if (modifiers == MODIF_ALT) {// this is reserved for manipulation
        SurfaceMeshPicker picker(camera());
        for(auto [name, mesh] : manager_SurfaceMesh){
            if(manager_bool[name+".pick"] == true){
                auto vert = picker.pick_vertex(mesh, x, y);
                if(vert.is_valid()){
                    unmark_select_point();
                    picking_mode = "faced_point";
                    mark_point(mesh->position(vert));
                    _picked_vert = vert; _picked_mesh = mesh;
                    auto pos = mesh->position(vert);
                    _mng_mani_0[name](_picked_vert.idx(),{pos.x, pos.y, pos.z}, _I_mat4_);
                    std::cout << mesh->position(vert) << std::endl;
                    std::cout << "Picking: " << vert.idx() << " " << name << std::endl; 
                    return false;
                }
            }
        }
        if(_picked_mesh != nullptr && button == BUTTON_RIGHT){
            unmark_select_point();
            return false;
        } 
    }
    if(modifiers == MODIF_CTRL){
        if(button == BUTTON_LEFT){
            _picked_polygon.push_back(vec2(static_cast<float>(x), static_cast<float>(y)));
            return false; // ? true or false
        }
        else if(button == BUTTON_RIGHT) {
            unmark_select_point();
        }
    };
    return Viewer::mouse_press_event(x, y, button, modifiers);
}

vector<vector<double>> worldMat(easy3d::mat4 mat, easy3d::vec3 rel){
    vector<vector<double>> ret(4, vector<double>(4, 0));
    for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++)
        ret[i][j] = mat(i, j);
    ret[0][3] = ret[0][3] - (mat(0, 0) * rel[0] + mat(0, 1) * rel[1] + mat(0, 2) * rel[2]);
    ret[1][3] = ret[1][3] - (mat(1, 0) * rel[0] + mat(1, 1) * rel[1] + mat(1, 2) * rel[2]);
    ret[2][3] = ret[2][3] - (mat(2, 0) * rel[0] + mat(2, 1) * rel[1] + mat(2, 2) * rel[2]);
    return ret;
}

template<typename R3>
R3 transform(R3 data, std::vector<std::vector<double>> trans, double weight) {
    assert(trans.size() == 4 && trans[0].size() == 4);
    double x = data[0];
    double y = data[1];
    double z = data[2];
    double x_ = trans[0][0] * x + trans[0][1] * y + trans[0][2] * z + trans[0][3] * weight;
    double y_ = trans[1][0] * x + trans[1][1] * y + trans[1][2] * z + trans[1][3] * weight;
    double z_ = trans[2][0] * x + trans[2][1] * y + trans[2][2] * z + trans[2][3] * weight;
    return { float(x_), float(y_), float(z_) };
}

bool ViewerImGui::mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) {
    if (modifiers == MODIF_SHIFT && _picked_pos) {
        ManipulatedFrame *frame = _picked_pos->manipulator()->frame();
        auto axis = ManipulatedFrame::NONE;
        if (pressed_key_ == KEY_X) axis = ManipulatedFrame::HORIZONTAL;
        else if (pressed_key_ == KEY_Y) axis = ManipulatedFrame::VERTICAL;
        else if (pressed_key_ == KEY_O) axis = ManipulatedFrame::ORTHOGONAL;
        switch (button) {
            case BUTTON_RIGHT:
                frame->action_rotate(x, y, dx, dy, camera_, axis);
                break;
            case BUTTON_LEFT:
                frame->action_translate(x, y, dx, dy, camera_, axis);
                break;
            default:
                break;
        }
        auto vert = _picked_pos->points()[0]; // Might be bug
        auto mat = worldMat(frame->worldMatrix(), vert);

        if(picking_mode == "faced_point"){
            auto pos = transform(vert, mat, 1); 
            auto name = _picked_mesh->name();
            if(_mng_mani_0.find(name) != _mng_mani_0.end())
                _mng_mani_0[name](_picked_vert.idx(),{pos.x, pos.y, pos.z}, mat);
        }
        else if(picking_mode == "point_cloud"){
            vector<vector<double>> _picked_points_pos_new;
            for(auto p : _picked_points_pos){
                auto pos_new = transform(p, mat, 1);
                _picked_points_pos_new.push_back(pos_new);
            }
            vec4 color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
            auto name = _picked_point_cloud_name;
            add_points("picked_point_clouds", _picked_points_pos_new, picked_point_size, color);
            update();
            if(_mng_mani_1.find(name) != _mng_mani_1.end())
                _mng_mani_1[name](_picked_points_idx, _picked_points_pos_new, mat);
        }
        return false;
    }
    if(modifiers == MODIF_CTRL){
        _picked_polygon.push_back(vec2(static_cast<float>(x), static_cast<float>(y)));
        return false;
    }
    else
        return Viewer::mouse_drag_event(x, y, dx, dy, button, modifiers);
}

int count_picked_clouds(PointCloud* cloud){
    int ret = 0;
    auto select = cloud->vertex_property<bool>("v:select");
    for(auto v : cloud->vertices())
        if(select[v]) ret++;
    return ret;
}

bool ViewerImGui::mouse_release_event(int x, int y, int button, int modifiers) {
    if (modifiers == MODIF_CTRL) {
        // _picked_polygon.push_back(_picked_polygon[0]);
        if (_picked_polygon.size() >= 3) {
            for(auto [name, cloud] : manager_PointCloud){
                if (cloud && manager_bool[name+".pick"] == true) {
                    PointCloudPicker picker(camera());
                    picker.pick_vertices(cloud, _picked_polygon, 0);
                    int n_picked = count_picked_clouds(cloud);
                    if(n_picked == 0) continue;
                    mark_point_clouds(cloud);
                    _picked_point_cloud_name = name;
                    std::cout << "Picked: " << n_picked << std::endl;
                }
            }
            _picked_polygon.clear();
        }
        return false;
    } 
    if(_picked_polygon.size() != 0) _picked_polygon.clear();
    return Viewer::mouse_release_event(x, y, button, modifiers);
}

void ViewerImGui::mark_point(easy3d::vec3 pos){
    vec4 color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
    add_mesh("picked_position", [&](int idx, int i){ return pos[i]; }, 1, 
                             [&](int idx, int i){return 0;}, 1, color, false);
    auto pts = manager_SurfaceMesh["picked_position"]->renderer()->points_drawables();
    for(auto p : pts) {
        p->set_point_size(picked_point_size);
        p->set_visible(true);
    };
    for(auto m : models()) if(m->name() == "picked_position")
        _picked_pos= m;
    if(!_picked_pos->manipulator()){
        _picked_pos->set_manipulator(new Manipulator(_picked_pos));
        _picked_pos->manipulator()->frame()->modified.connect(this,
                static_cast<void (ViewerImGui::*)(void) const>(&ViewerImGui::update));
    }
    update();
}

int ViewerImGui::mark_point_clouds(PointCloud* cloud){
    unmark_select_point();
    picking_mode = "point_cloud";
    vec4 color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
    vector<vector<double>> pts;
    vec3 center = {0, 0, 0};
    auto select = cloud->vertex_property<bool>("v:select");
    _picked_points_idx.clear();
    _picked_points_pos.clear();
    for(auto v : cloud->vertices()){
        auto pos = cloud->position(v);
        if(select[v]){
            cloud->vertex_property<bool>("v:select")[v] = false;
            center += pos;
            _picked_points_idx.push_back(v.idx());
            _picked_points_pos.push_back({pos[0], pos[1], pos[2]});
        }
    }
    if(_picked_points_idx.size() == 0){
        unmark_select_point(); return 0;
    }
    center /= _picked_points_idx.size();
    add_points("picked_point_clouds", _picked_points_pos, picked_point_size, color);
    mark_point(center);
    _picked_pos->renderer()->set_visible(false);
    update();
    return _picked_points_idx.size();
}

void ViewerImGui::unmark_select_point(){
    if(picking_mode == "faced_point"){
        add_mesh("picked_position", {}, 0, {}, 0, {}, false);
        _picked_pos = NULL;
    }
    else if(picking_mode == "point_cloud"){
        vector<vector<double>> no_pts;
        _picked_points_idx.clear();
        _picked_points_pos.clear();
        add_points("picked_point_clouds", no_pts, picked_point_size, {});
        add_mesh("picked_position", {}, 0, {}, 0, {}, false);
        _picked_pos = NULL;
    }
    picking_mode = "none";
    update();
}

void ViewerImGui::draw() const {
    // default draw

    // transparency draw;
    std::vector<PointsDrawable*> solid_points;
    std::vector<LinesDrawable*> solid_lines;
    std::vector<TrianglesDrawable*> solid_faces;
    std::vector<TrianglesDrawable*> transparent_faces;
    for(auto m : models()) {
        if(!m->renderer()->is_visible()) continue;
        auto name_transparency = m->name() + ".transparency";
        auto name_line_vis = m->name() + ".edge_visible";
        if(manager_bool[name_transparency]) {
            for(auto f : m->renderer()->triangles_drawables()) if(f->is_visible())
                transparent_faces.push_back(dynamic_cast<TrianglesDrawable*>(f));
        } 
        else {
            for(auto f : m->renderer()->triangles_drawables()) if(f->is_visible())
                solid_faces.push_back(dynamic_cast<TrianglesDrawable*>(f));
        }
        for(auto l : m->renderer()->lines_drawables()) if(l->is_visible() && manager_bool[name_line_vis])
            solid_lines.push_back(dynamic_cast<LinesDrawable*>(l));
        for(auto p : m->renderer()->points_drawables()) if(p->is_visible()){
            solid_points.push_back(dynamic_cast<PointsDrawable*>(p));
        }
    }

    __transparancy_render__->draw(transparent_faces);
    for(auto l : solid_lines) l->draw(camera());
    for(auto p : solid_points) p->draw(camera());

    for(auto [name, mat] : manager_axes_wordMat){
        Camera c(*camera()); 
        c.set_modelview_matrix(camera()->modelViewMatrix() * mat);
        drawable_axes_->draw(&c);
    }

    if(solid_points.size() != 0 || solid_lines.size() != 0){
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(0.5f, -0.0001f);
        for(auto f : solid_faces) f->draw(camera());
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    else {
        for(auto f : solid_faces) f->draw(camera());
    }
    for (auto d : drawables_) {
        if (d->is_visible())
            d->draw(camera());
    }

    // Viewer::draw();
}

// draw imgui and so on
void ViewerImGui::post_draw() {
    custom_imgui();
    static bool show_about = false;
    if (show_about) {
        ImGui::SetNextWindowPos(ImVec2(static_cast<float>(width()) * 0.5f, static_cast<float>(height()) * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::Begin("About Easy3D ImGui Viewer", &show_about, ImGuiWindowFlags_NoResize);
        ImGui::Text("This viewer shows how to use ImGui for GUI creation and event handling");
        ImGui::Separator();
        ImGui::Text(
            "\n"
            "Liangliang Nan\n"
            "liangliang.nan@gmail.com\n"
            "https://3d.bk.tudelft.nl/liangliang/\n"
        );
        ImGui::End();
    }
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5, 8));
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("Open", "Ctrl+O"))
                open();
            if (ImGui::MenuItem("Save As...", "Ctrl+S"))
                save();
            ImGui::Separator();
            if (ImGui::MenuItem("Quit", "Alt+F4"))
                glfwSetWindowShouldClose(window_, GLFW_TRUE);
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("View"))
        {
            if (ImGui::MenuItem("Snapshot", nullptr))
                snapshot(); // TODO scaling
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Help")) {
            ImGui::MenuItem("About", nullptr, &show_about);
            ImGui::EndMenu();
        }
        menu_height_ = ImGui::GetWindowHeight();
        ImGui::EndMainMenuBar();
    }
    ImGui::PopStyleVar();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); 
    auto texter = texter_;

    texter_ = nullptr;
    if (_picked_polygon.size() >= 3) {
        // draw the boundary of the rect/lasso
        shape::draw_polygon_wire(_picked_polygon, vec4(1.0f, 0.0f, 0.0f, 1.0f), width(), height(), -1.0f);
        // draw its transparent face
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        shape::draw_polygon_filled(_picked_polygon, vec4(1.0f, 0.0f, 0.0f, 0.2f), width(), height(), -0.9f);
        glDisable(GL_BLEND);
    }

    // Viewer::post_draw();
    texter_ = texter;
    draw_corner_axes();

}

void ViewerImGui::set_manager_bool(std::string name, bool v){
    manager_bool[name] = v;
    // update manager_bool with fact;
    auto split = [&](std::string s, char t){
        std::vector<std::string> ret({{}});
        for(auto c : s){
            if(c != t) ret.back().push_back(c);
            else ret.push_back({});
        }
        return ret;
    };
    auto names = split(name, '.');
    for(auto m : models()) if(m->name() == names[0])
    if(names[1] == "edge_visible")
        m->renderer()->get_lines_drawable("edges")->set_visible(v);
}
}
