from tkinter import *
import matplotlib.pyplot as plt
import model

window = Tk()

TimeStart = DoubleVar()
TimeEnd = DoubleVar()
TW = DoubleVar()
Rk = DoubleVar()
Le = DoubleVar()
Rad = DoubleVar()
Lk = DoubleVar()
Ck = DoubleVar()
Uc0 = DoubleVar()
I0 = DoubleVar()
Tau = DoubleVar()
M = IntVar()
UseLamp = BooleanVar()


def build():
    model.Mode = M.get()
    model.TimeStep = Tau.get()
    model.CurrentDefault = I0.get()
    model.VoltageDefault = Uc0.get()
    model.CapacityDefault = Ck.get()
    model.InductanceDefault = Lk.get()
    model.LampRadiusDefault = Rad.get()
    model.LengthBetweenElectrodesDefault = Le.get()
    model.ResistanceDefault = Rk.get()
    model.EdgeTemperature = TW.get()
    model.UsingGas = UseLamp.get()
    model.TimeDefault = TimeStart.get()
    model.TimeEndDefault = TimeEnd.get()
    if model.Mode == 1:
        r = model.runge_kutta_4()
        show(r)
    elif model.Mode == 2:
        r = model.runge_kutta_2()
        show(r)
    else:
        r = model.runge_kutta_2(alpha=0.5)
        show(r)


def configure_mode(root):
    frame_mode = LabelFrame(root, text="Запуск")
    frame_mode.grid(row=7, column=4, rowspan=5, columnspan=3)
    M.set(model.Mode)
    UseLamp.set(model.UsingGas)
    mode_txt = Radiobutton(frame_mode, text="RK4", variable=M, value=1)
    mode_txt.grid(row=1, column=0)

    moderk2_txt = Radiobutton(frame_mode, text="RK2", variable=M, value=2)
    moderk2_txt.grid(row=1, column=1)

    moderk2mid_txt = Radiobutton(frame_mode, text="RK2mid", variable=M, value=3)
    moderk2mid_txt.grid(row=1, column=2)

    use = Checkbutton(frame_mode, text="С трубкой?", variable=UseLamp, onvalue=True, offvalue=False)
    use.grid(row=2, column=0)
    show_button = Button(frame_mode, text="Построить графики", command=build)
    show_button.grid(row=3, column=0)


def configure_time(root):
    frame_start = LabelFrame(root, text="Начальные условия")
    frame_start.grid(row=12, column=0, rowspan=5, columnspan=3)

    TimeStart.set(model.TimeDefault)
    TimeEnd.set(model.TimeEndDefault)
    u0_lbl = Label(frame_start, text="Tbeg")
    u0_lbl.grid(row=2, column=0)
    u0_txt = Entry(frame_start, textvariable=TimeStart)
    u0_txt.grid(row=2, column=1)
    i0_lbl = Label(frame_start, text="Tend")
    i0_lbl.grid(row=3, column=0)
    i0_txt = Entry(frame_start, textvariable=TimeEnd)
    i0_txt.grid(row=3, column=1)


def configure_start(root):
    frame_start = LabelFrame(root, text="Начальные условия")
    frame_start.grid(row=7, column=0, rowspan=5, columnspan=3)
    Tau.set(model.TimeStep)
    Uc0.set(model.VoltageDefault)
    I0.set(model.CurrentDefault)
    tau_lbl = Label(frame_start, text="Tau")
    tau_lbl.grid(row=1, column=0)
    tau_txt = Entry(frame_start, textvariable=Tau)
    tau_txt.grid(row=1, column=1)
    u0_lbl = Label(frame_start, text="U0")
    u0_lbl.grid(row=2, column=0)
    u0_txt = Entry(frame_start, textvariable=Uc0)
    u0_txt.grid(row=2, column=1)
    i0_lbl = Label(frame_start, text="I0")
    i0_lbl.grid(row=3, column=0)
    i0_txt = Entry(frame_start, textvariable=I0)
    i0_txt.grid(row=3, column=1)


def configure_freq(root):
    frame_freq = LabelFrame(root, text="Колебательный контур")
    frame_freq.grid(row=4, column=4, rowspan=3, columnspan=3)
    Lk.set(model.InductanceDefault)
    Ck.set(model.CapacityDefault)
    lk_lbl = Label(frame_freq, text="Lk")
    lk_lbl.grid(row=1, column=0)
    ck_lbl = Label(frame_freq, text="Ck")
    ck_lbl.grid(row=2, column=0)
    lk_txt = Entry(frame_freq, textvariable=Lk)
    lk_txt.grid(row=1, column=1)
    ck_txt = Entry(frame_freq, textvariable=Ck)
    ck_txt.grid(row=2, column=1)


def configure_size(root):
    frame_size = LabelFrame(root, text="Размеры трубки")
    frame_size.grid(row=4, column=0, rowspan=3, columnspan=3)
    Le.set(model.LengthBetweenElectrodesDefault)
    Rad.set(model.LampRadiusDefault)
    le_lbl = Label(frame_size, text="L")
    le_lbl.grid(row=2, column=0)
    le_text = Entry(frame_size, textvariable=Le)
    le_text.grid(row=2, column=1)
    rad_lbl = Label(frame_size, text="R")
    rad_lbl.grid(row=1, column=0)
    rad_text = Entry(frame_size, textvariable=Rad)
    rad_text.grid(row=1, column=1)


def configure_resistance(root):
    frame_resistance = LabelFrame(root, text="Постоянное сопротивление")
    frame_resistance.grid(row=0, column=4, rowspan=3, columnspan=3)

    Rk.set(model.ResistanceDefault)
    rk_text = Entry(frame_resistance, textvariable=Rk)
    rk_text.grid(row=1, column=1)


def configure_temperature_gui(root):
    frame_temperature = LabelFrame(root, text="Краевая температура")
    frame_temperature.grid(row=0, column=0, rowspan=3, columnspan=3)
    TW.set(model.EdgeTemperature)
    temperature_end_text = Entry(frame_temperature, textvariable=TW)
    temperature_end_text.grid(row=1, column=1)


def configure_gui(root):
    root.title("Lab 2")
    configure_temperature_gui(root)
    configure_resistance(root)
    configure_size(root)
    configure_freq(root)
    configure_start(root)
    configure_mode(root)
    configure_time(root)


def show(result):
    plt.title('mode = ' + str(model.Mode))
    plt.subplot(3, 3, 1)
    plt.plot(result[0], result[1])
    plt.ylabel("I, А")

    plt.subplot(3, 3, 4)
    plt.plot(result[0], result[2])
    plt.ylabel("Uc, В")

    plt.subplot(3, 3, 7)
    plt.plot(result[0], result[3])
    plt.ylabel("Rp, Ом")

    plt.subplot(3, 3, 3)
    plt.plot(result[0], result[4])
    plt.ylabel("Т0(0), К")

    us = []
    for i in range(len(result[1])):
        us.append(result[1][i] * result[3][i])

    plt.subplot(3, 3, 6)
    plt.plot(result[0], us)
    plt.ylabel("I0*Rp, В")
    plt.xlabel("t, c")
    plt.show()


configure_gui(window)
window.mainloop()
