"""NASA Python unit helpers; scalar and NumPy inputs keep their shape."""


def celsius_to_kelvin(value):
    return value + 273.15


def kelvin_to_celsius(value):
    return value - 273.15


def fahrenheit_to_kelvin(value):
    return (value - 32) / 1.8 + 273.15


def kelvin_to_fahrenheit(value):
    return (value - 273.15) * 1.8 + 32


def rankine_to_kelvin(value):
    return value / 1.8


def kelvin_to_rankine(value):
    return value * 1.8


def atm_to_bar(value):
    return value * 1.01325


def bar_to_atm(value):
    return value / 1.01325


def psi_to_bar(value):
    return value * (1.01325 / 14.695948775513)


def bar_to_psi(value):
    return value / (1.01325 / 14.695948775513)


def mmhg_to_bar(value):
    return value * (1.01325 / 760)


def bar_to_mmhg(value):
    return value / (1.01325 / 760)


def g_per_cm3_to_kg_per_m3(value):
    return value * 1000


def kg_per_m3_to_g_per_cm3(value):
    return value / 1000


def cm3_per_g_to_m3_per_kg(value):
    return value / 1000


def m3_per_kg_to_cm3_per_g(value):
    return value * 1000


def m_per_s_to_ft_per_s(value):
    return value * 3.28084


def ft_per_s_to_m_per_s(value):
    return value / 3.28084


def joule_to_cal(value):
    return value / 4.184


def cal_to_joule(value):
    return value * 4.184


def kj_per_kg_to_cal_per_g(value):
    return value / 4.184


def cal_per_g_to_kj_per_kg(value):
    return value * 4.184


__all__ = [name for name in globals() if not name.startswith("_")]
