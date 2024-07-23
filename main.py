import numpy as np
import dearpygui.dearpygui as dpg

global selected_point, is_drag
is_drag = False
selected_point = None


dpg.create_context()
DEG = 3
WIDTH = 500
HEIGHT = 500
R = 70
global cps, method, t_selected
cps = np.array([[100, 400], [200, 300], [300, 300], [400, 400]])
t_selected = 0.5
method = "Polynomial"
methods = ["Polynomial", "Bernstein Poly", "de Casteljau", "Seiler (diff. terms)", "Seiler (pure lerp)", "Seiler (offsets)"]
dpg.add_window(no_scrollbar=True, tag="main_window")
dpg.add_window(no_scrollbar=True, tag="settings")

def set_method(sender, app_data):
    global method
    method = app_data
    draw()

def set_t(sender, app_data):
    global t_selected
    t_selected = app_data
    draw()

dpg.add_slider_float(label="t", min_value=0.0, max_value=1.0, callback=set_t, parent="settings", default_value=t_selected)
dpg.add_text("methods", parent="settings")
dpg.add_radio_button(methods, callback=set_method, parent="settings")

def factorial(n):
    result = 1
    for i in range(n):
        result *= i+1
    return result

def binomial_coefficients(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def de_casteljau(t, bs):
    if len(bs) == 1:
        return bs[0]
    result = []
    for i, b in enumerate(bs):
        if i == 0:
            continue
        result.append(L(bs[i-1], b, t))
    result = de_casteljau(t, result)
    return result

# def de_casteljau(t, coefs):
#     beta = coefs.copy()
#     n = len(beta)
#     for j in range(1, n):
#         for k in range(n - j):
#             beta[k] = beta[k] * (1 - t) + beta[k + 1] * t
#     return beta[0]

def L(a, b, t):
    return (1 - t) * a + t * b

def get_ds():
    ds = {}
    ds[0] = cps[0]
    ds[DEG] = cps[-1]
    ds[1] = DEG * (cps[1] - cps[0]) - (cps[-1] - cps[0])
    if ds.get(DEG - 1) is None:
        ds[DEG - 1] = DEG * (cps[-2] - cps[-1]) - (cps[0] - cps[-1])
    if ds.get(2) is None:
        ds[2] = binomial_coefficients(DEG, 2) * (cps[2] - cps[1]) - binomial_coefficients(DEG - 2, 2) * (cps[1] - cps[0]) - (DEG - 3) * (cps[-2] - cps[-1]) - 3 * (cps[-2] - cps[1])
    if ds.get(DEG - 2) is None:
        ds[DEG - 2] = binomial_coefficients(DEG, 2) * (cps[-3] - cps[-2]) - binomial_coefficients(DEG - 2, 2) * (cps[-2] - cps[-1]) - (DEG - 3) * (cps[1] - cps[0]) - 3 * (cps[1] - cps[-2])
    return ds

# coeff of Polynomial form. Ref: https://en.wikipedia.org/wiki/B%C3%A9zier_curve
def C(j):
    k = 1
    for m in range(j):
        k *= (DEG - m)
    s = sum([(-1) ** (i + j) * cps[i] / (factorial(i) * factorial(j-i)) for i in range(j+1)])
    return k * s


def D(i, t, ds):
    if i * 2 == DEG + 1:
        return 0
    if i * 2 == DEG:
        return ds[i]
    return L(ds[i], ds[DEG - i], t) + (1 - t) * t * D(i+1, t, ds)

def E(i, t, ds, positive=True):
    if DEG - 1 <= 2 * i and 2 * i <= DEG + 1:
        return ds[i]
    i_next = i + 1 if positive else i - 1
    return ds[i] + (1 - t) * t * E(i_next, t, ds, positive)

gray = [120, 120, 120]
yellow = [200, 200, 0]
blue = [120, 120, 255]
red = [255, 0, 0]
green = [0, 255, 0]

# in case of degree = 3
def draw_cubic():
    # init canvas
    if dpg.does_item_exist("canvas"):
        dpg.delete_item("canvas")
    dpg.add_drawlist(width=WIDTH, height=HEIGHT, tag="canvas", parent="main_window")

    # draw control points
    for i, cp in enumerate(cps):
        dpg.draw_circle(cp, 5, parent="canvas", color=gray)
        if i > 0:
            dpg.draw_line(cp, cps[i-1], thickness=1, parent="canvas", color=gray)

    # draw curve
    N = 100
    if method == "Polynomial":
        p0 = cps[0]
        p1 = 3 * (cps[1] - cps[0])
        p2 = 3 * (cps[0] - 2 * cps[1] + cps[2])
        p3 = -cps[0] + 3 * cps[1] - 3 * cps[2] + cps[3]
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            current = t ** 3 * p3 + t ** 2 * p2 + t * p1 + p0
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    elif method == "Bernstein Poly":
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            current = (1-t) ** 3 * cps[0] + 3 * (1-t) ** 2 * t * cps[1] + 3 * (1-t) * t ** 2 * cps[2] + t ** 3 * cps[3]
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    elif method == "de Casteljau":
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            b01 = L(cps[0], cps[1], t)
            b12 = L(cps[1], cps[2], t)
            b23 = L(cps[2], cps[3], t)
            b02 = L(b01, b12, t)
            b13 = L(b12, b23, t)
            current = L(b02, b13, t)
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    else:
        d1 = 3 * (cps[1] - cps[0]) - (cps[3] - cps[0])
        d2 = 3 * (cps[2] - cps[3]) - (cps[0] - cps[3])
        if method == "Seiler (diff. terms)":
            for i in range(N+1):
                t = i / N
                if i == 0:
                    prev = cps[0]
                b03 = L(cps[0], cps[3], t)
                d12 = L(d1, d2, t)
                current = b03 + (1 - t) * t * d12
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current
        elif method == "Seiler (pure lerp)":
            s1 = d1 + cps[0]
            s2 = d2 + cps[3]
            for i in range(N+1):
                t = i / N
                if i == 0:
                    prev = cps[0]
                s12 = L(s1, s2, t)
                b03 = L(cps[0], cps[3], t)
                current = L(b03, s12, (1 - t) * t)
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current
        elif method == "Seiler (offsets)":
            for i in range(N+1):
                if i == 0:
                    prev = cps[0]
                t = i / N
                e1 = cps[0] + (1 - t) * t * d1
                e2 = cps[3] + (1 - t) * t * d2
                current = L(e1, e2, t)
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current
    
    # draw auxiliary lines and circles
    t = t_selected
    if method == "Polynomial":
        p0 = cps[0]
        p1 = 3 * (cps[1] - cps[0])
        p2 = 3 * (cps[0] - 2 * cps[1] + cps[2])
        p3 = -cps[0] + 3 * cps[1] - 3 * cps[2] + cps[3]
        current = t ** 3 * p3 + t ** 2 * p2 + t * p1 + p0
    elif method == "Bernstein Poly":
        current = (1-t) ** 3 * cps[0] + 3 * (1-t) ** 2 * t * cps[1] + 3 * (1-t) * t ** 2 * cps[2] + t ** 3 * cps[3]
    elif method == "de Casteljau":
        b01 = L(cps[0], cps[1], t)
        b12 = L(cps[1], cps[2], t)
        b23 = L(cps[2], cps[3], t)
        b02 = L(b01, b12, t)
        b13 = L(b12, b23, t)
        dpg.draw_line(b01, b12, thickness=1, parent="canvas", color=yellow)
        dpg.draw_line(b12, b23, thickness=1, parent="canvas", color=yellow)
        dpg.draw_line(b13, b02, thickness=1, parent="canvas", color=yellow)
        dpg.draw_circle(b01, 5, parent="canvas", color=yellow)
        dpg.draw_circle(b12, 5, parent="canvas", color=yellow)
        dpg.draw_circle(b23, 5, parent="canvas", color=yellow)
        dpg.draw_circle(b02, 5, parent="canvas", color=yellow)
        dpg.draw_circle(b13, 5, parent="canvas", color=yellow)
        current = L(b02, b13, t)
    else:
        d1 = 3 * (cps[1] - cps[0]) - (cps[3] - cps[0])
        d2 = 3 * (cps[2] - cps[3]) - (cps[0] - cps[3])
        b03 = L(cps[0], cps[3], t)
        s1 = d1 + cps[0]
        s2 = d2 + cps[3]
        if method == "Seiler (diff. terms)":
            dpg.draw_line(cps[0], cps[3], thickness=1, parent="canvas", color=blue)
            d12 = L(d1, d2, t)
            current = b03 + (1 - t) * t * d12
            dpg.draw_circle(b03, 5, parent="canvas", color=blue)
            dpg.draw_line(b03, d12 + b03, thickness=1, parent="canvas", color=blue)
        elif method == "Seiler (pure lerp)":
            dpg.draw_line(cps[0], cps[3], thickness=1, parent="canvas", color=blue)
            dpg.draw_line(s1, cps[0], thickness=1, parent="canvas", color=green)
            dpg.draw_line(s2, cps[3], thickness=1, parent="canvas", color=green)
            dpg.draw_line(s2, s1, thickness=1, parent="canvas", color=green)
            s12 = L(s1, s2, t)
            current = L(b03, s12, (1 - t) * t)
            dpg.draw_circle(b03, 5, parent="canvas", color=blue)
            dpg.draw_circle(s12, 5, parent="canvas", color=green)
            dpg.draw_line(b03, s12, thickness=1, parent="canvas", color=blue)
        elif method == "Seiler (offsets)":
            dpg.draw_line(s1, cps[0], thickness=1, parent="canvas", color=green)
            dpg.draw_line(s2, cps[3], thickness=1, parent="canvas", color=green)
            e1 = cps[0] + (1 - t) * t * d1
            e2 = cps[3] + (1 - t) * t * d2
            dpg.draw_circle(e1, 5, parent="canvas", color=red)
            dpg.draw_circle(e2, 5, parent="canvas", color=red)
            dpg.draw_line(e1, e2, thickness=1, parent="canvas", color=red)
            current = L(e1, e2, t)
    dpg.draw_circle(current, 5, parent="canvas")

# in case of 2 <= degree <= 5
def draw():
    # init canvas
    if dpg.does_item_exist("canvas"):
        dpg.delete_item("canvas")
    dpg.add_drawlist(width=WIDTH, height=HEIGHT, tag="canvas", parent="main_window")

    # draw control points
    for i, cp in enumerate(cps):
        dpg.draw_circle(cp, 5, parent="canvas", color=gray)
        if i > 0:
            dpg.draw_line(cp, cps[i-1], thickness=1, parent="canvas", color=gray)

    # draw curve
    N = 100
    if method == "Polynomial":
        coeffs = []
        for i in range(DEG+1):
            coeffs.append(C(i))
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            current = 0
            for j in range(DEG+1):
                current += t ** j * coeffs[j]
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    elif method == "Bernstein Poly":
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            current = 0
            for j, cp in enumerate(cps):
                coeff = binomial_coefficients(DEG, j)
                current += coeff * t ** j * (1 - t) ** (DEG - j) * cp
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    elif method == "de Casteljau":
        for i in range(N+1):
            t = i / N
            if i == 0:
                prev = cps[0]
            current = de_casteljau(t, cps)
            dpg.draw_line(prev, current, thickness=1, parent="canvas")
            prev = current
    else:
        ds = get_ds()
        if method == "Seiler (diff. terms)":
            for i in range(N+1):
                t = i / N
                if i == 0:
                    prev = cps[0]
                current = L(cps[0], cps[-1], t) + (1 - t) * t * D(1, t, ds)
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current
        elif method == "Seiler (pure lerp)":
            s1 = ds[1] + cps[0]
            s2 = ds[DEG - 1] + cps[3]
            for i in range(N+1):
                t = i / N
                if i == 0:
                    prev = cps[0]
                s12 = L(s1, s2, t)
                b03 = L(cps[0], cps[3], t)
                current = L(b03, s12, (1 - t) * t)
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current
        elif method == "Seiler (offsets)":
            ds = get_ds()
            for i in range(N+1):
                if i == 0:
                    prev = cps[0]
                t = i / N
                current = L(E(0, t, ds), E(DEG, t, ds, False), t)
                dpg.draw_line(prev, current, thickness=1, parent="canvas")
                prev = current

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

draw()
dpg.create_viewport()
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()