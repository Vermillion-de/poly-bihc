#pragma once
#include <easy3d/viewer/viewer.h>
#include <easy3d/util/initializer.h>
#include <easy3d/util/dialog.h>
#include <easy3d/util/file_system.h>
#include <functional>

struct ImGuiContext;

namespace easy3d {
using std::string;
using std::vector;
using std::function;
const vec4 _green{1.0f, 0.8f, 0.4f, 1.0f};
const vec4 _black{0.0f, 0.0f, 0.0f, 1.0f};
	class Model;
	class PointCloud;
	void init();
    class ViewerImGui : public Viewer
	{
	public:
        explicit ViewerImGui(
            const string& title = "Viewer GUI", int samples = 4,
            int gl_major = 3,  int gl_minor = 2,  bool full_screen = false, bool resizable = true,
			int depth_bits = 24, int stencil_bits = 8, int width = 2400, int height = 1600 
		);
        ~ViewerImGui() override;

    public: 
		// window manager of parameters & surface
        // void manage_float(string, float*, function<void()> call_back = [](){});
        void manage_float(string, float*, float f_min=0, float f_max=1, function<void()> call_back = [](){});
        void manage_int(string, int*, int i_min=0, int i_max=10, function<void()> call_back = [](){});
		void manage_int2(string ,int*, int i_min=0, int i_max=10, function<void()> call_back = [](){});
        void manage_clear(string);
		void manage_clear(string, string, string);
        void manage_button(string, function<void()>);
        void manage_checkbox(string, bool*, function<void()> callback=[](){});
        void manage_keycallback(char, function<void()>);
        void manage_selectbox(string, vector<string>, string*, function<void()> call_back = [](){});
        void manage_selectradio(string, vector<string>, string*, function<void()> call_back = [](){});

		void manage_picked_faced_point(string, function<void(int, vector<double>, vector<vector<double>> rot)>);
		void manage_picked_point_cloud(string, function<void(vector<int>, vector<vector<double>>, vector<vector<double>> rot)>);

		void set_manager_bool(string, bool);

		// add geometric permitives
        void add_handle_surf(string name, function<vector<double>(double, double)> handle, vector<double> u_range, vector<double> v_range, int res=100, vec4 color=_green);
		void add_points(string name, vector<vector<double>> pts, float point_size=50, vec4 color=_green);
		void add_mesh(string name, function<float(int, int)> f_vert, int nv, function<int(int, int)> f_face, int nf, vec4 color=_green, bool smooth_shading=true);
		void add_graph(string name, vector<vector<double>> pts, vector<vector<int>> edges, float line_width=10, float point_size = 20, vec4 pts_color=_green, vec4 edge_color=_black);
		void add_implot_lines(string name, vector<vector<vector<double>>> lines, float line_width=10);
		void set_axes(string name, vector<vector<double>> word_mat);

        void custom_imgui();
        void render_settings();

		void log_imgui(string log);

		void save_camera_parameter(string path);

		// count for object manipulation
		bool mouse_press_event(int x, int y, int button, int modifiers) override;
		bool mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) override;
		bool mouse_release_event(int x, int y, int button, int modifiers) override;
		void mark_point(vec3 pos);
		int mark_point_clouds(PointCloud*);
		void unmark_select_point();
		void draw() const override;

	protected: // No changed lifetime
		void init() override; 
		void reload_font(int font_size = 16);
		void pre_draw() override; // draw the widgets
		void post_draw() override; //  the widgets
		void post_resize(int, int) override; //  the widgets

        // callbacks
		bool callback_event_cursor_pos(double x, double y) override;
		bool callback_event_mouse_button(int button, int action, int modifiers) override;
		bool callback_event_keyboard(int key, int action, int modifiers) override;
		bool callback_event_character(unsigned int codepoint) override;
		bool callback_event_scroll(double dx, double dy) override;

	public:
		float picked_point_size = 30;

	protected:
		float refresh_rate = 60;
        static ImGuiContext *	context_;// Single global context by default, but can be overridden by the user
		float menu_height_;
	};
}