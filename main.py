import numpy as np
import dearpygui.dearpygui as dpg

global selected_point, is_drag
is_drag = False
selected_point = None


def binomial_coefficients(n, k):
    a, b, c = 1, 1, 1
    for i in range(n):
        a *= i+1
    for i in range(k):
        b *= i+1
    for i in range(n-k):
        c *= i+1
    return a / (b * c)

def L(a, b, t):
    return (1 - t) * a + t * b

def update():
    N = 100
    points = []
    if method == "Polynomial":
        p0 = cps[0]
        p1 = 3 * (cps[1] - cps[0])
        p2 = 3 * (cps[0] - 2 * cps[1] + cps[2])
        p3 = -cps[0] + 3 * cps[1] - 3 * cps[2] + cps[3]
        for i in range(N+1):
            t = i / N
            p = t ** 3 * p3 + t ** 2 * p2 + t * p1 + p0
            points.append(p)
    elif method == "Bernstein Poly":
        for i in range(N+1):
            t = i / N
            # current = 0
            # for j, cp in enumerate(cps):
            #     coeff = binomial_coefficients(DEG, j)
            #     current += coeff * t ** j * (1 - t) ** (DEG - j) * cp
            p = (1-t) ** 3 * cps[0] + 3 * (1-t) ** 2 * t * cps[1] + 3 * (1-t) * t ** 2 * cps[2] + t ** 3 * cps[3]
            points.append(p)
    elif method == "de Casteljau":
        for i in range(N+1):
            t = i / N
            b01 = L(cps[0], cps[1], t)
            b12 = L(cps[1], cps[2], t)
            b23 = L(cps[2], cps[3], t)
            b02 = L(b01, b12, t)
            b13 = L(b12, b23, t)
            p = L(b02, b13, t)
            points.append(p)
    elif method == "Seiler (diff. terms)":
        d1 = 3 * (cps[1] - cps[0]) - (cps[3] - cps[0])
        d2 = 3 * (cps[2] - cps[3]) - (cps[0] - cps[3])
        for i in range(N+1):
            t = i / N
            b03 = L(cps[0], cps[3], t)
            d12 = L(d1, d2, t)
            p = b03 + (1 - t) * t * d12
            points.append(p)
    elif method == "Seiler (pure lerp)":
        d1 = 3 * (cps[1] - cps[0]) - (cps[3] - cps[0])
        d2 = 3 * (cps[2] - cps[3]) - (cps[0] - cps[3])
        s1 = d1 + cps[0]
        s2 = d2 + cps[3]
        for i in range(N+1):
            t = i / N
            s12 = L(s1, s2, t)
            b03 = L(cps[0], cps[3], t)
            p = L(b03, s12, (1 - t) * t)
            points.append(p)
    elif method == "Seiler (offsets)":
        d1 = 3 * (cps[1] - cps[0]) - (cps[3] - cps[0])
        d2 = 3 * (cps[2] - cps[3]) - (cps[0] - cps[3])
        for i in range(N+1):
            t = i / N
            e1 = cps[0] + (1 - t) * t * d1
            e2 = cps[3] + (1 - t) * t * d2
            p = L(e1, e2, t)
            points.append(p)
    return points

def draw():
    if dpg.does_item_exist("canvas"):
        dpg.delete_item("canvas")
    dpg.add_drawlist(width=WIDTH, height=HEIGHT, tag="canvas", parent="main_window")
    _ = [dpg.draw_circle(cp, 5, parent="canvas") for cp in cps]

    points = update()
    for i, current in enumerate(points):
        if i == 0:
            prev = cps[0]
        dpg.draw_line(prev, current, thickness=1, parent="canvas")
        prev = current

def set_method(sender, app_data):
    global method
    method = app_data
    draw()

dpg.create_context()
DEG = 3
WIDTH = 500
HEIGHT = 500
R = 70
global cps, method
cps = np.array([[100, 100], [200, 200], [300, 300], [400, 400]])
method = "Polynomial"
methods = ["Polynomial", "Bernstein Poly", "de Casteljau", "Seiler (diff. terms)", "Seiler (pure lerp)", "Seiler (offsets)"]
dpg.add_window(no_scrollbar=True, tag="main_window")
dpg.add_window(no_scrollbar=True, tag="settings")
dpg.add_radio_button(methods, callback=set_method, parent="settings")
draw()

def mouse_move(sender, app_data):
    if dpg.get_item_alias(dpg.get_active_window()) == "main_window":
        global selected_point
        if is_drag:
            if selected_point is not None:
                global cps
                pos = np.array(dpg.get_mouse_pos(local=True))
                cps[selected_point] = pos
                draw()
        else:
            pos = np.array(dpg.get_mouse_pos(local=True))
            dists = [np.linalg.norm(pos-p) for p in cps]
            nearest_point = dists.index(min(dists))
            if dists[nearest_point] > 50:
                selected_point = None
                dpg.delete_item("selected_point")
            elif nearest_point != selected_point:
                selected_point = nearest_point
                dpg.delete_item("selected_point")
                dpg.draw_circle(cps[selected_point], 5, color=(255, 0, 0), parent="canvas", tag="selected_point")

def mouse_click(sender, app_data):
    global is_drag
    is_drag = True

def mouse_release(sender, app_data):
    global is_drag
    is_drag = False

with dpg.handler_registry():
    dpg.add_mouse_click_handler(button=dpg.mvMouseButton_Left, callback=mouse_click)
    dpg.add_mouse_release_handler(button=dpg.mvMouseButton_Left, callback=mouse_release)
    dpg.add_mouse_move_handler(callback=mouse_move)

dpg.create_viewport()
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()