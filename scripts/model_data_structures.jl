
"""
QData
structure to hold the flow rate data and corresponding end times.
"""
struct QData{QT, T}
    Q::QT # flow velocity in m/s
    end_times::T
end

"""
    VData{QT, T}
Data structure to hold the flow velocity data and corresponding end times.
"""
struct VData{QT, T}
    v::QT # flow velocity in m/s
    end_times::T
end

"""
    VData{QT, T}
Data structure to hold the flow velocity data and corresponding start times.
"""
struct VDataS{QT, T}
    v::QT # flow velocity in m/s
    start_times::T
end

"""
    VDataA{QT, T}
Data structure to hold the flow velocity data and corresponding start and end times.
"""
struct VDataA{QT, T}
    v::QT # flow velocity in m/s
    start_times::T
    end_times::T
end

"""
    CinData{QT, T}
Data structure to hold inflow concentration and time.
This structure holds the inflow concentration data and the corresponding time.
"""
struct CinData{QT, T}
    c_in::QT # inflow concentration in mM
    t_in::T # time of the inflow concentration
end

"""
    conc_ds{C, T}
Data structure to hold concentration data and time.
This structure holds the concentration data and the corresponding time.
"""
struct conc_ds{C, T}
    conc::C # concentration in mM
    t::T # time in seconds since t0
    conc_ds(conc::C, t::T) where {C,T} = size(conc) == size(t) ? new{C,T}(conc, t) : error("Size mismatch between concentration and time")
end

"""
    ds{I, D}
Data structure to hold the dataset for each column.
This structure holds the column index and the concentration datasets for:
    NO3, NO2, SO4, Fe, pH, and EC.
"""
struct ds{I, D}
    column::I
    no3::D
    no2::D
    so4::D
    fe::D
    pH::D
    ec::D
end

"""
    ds{I, D}
Data structure to hold the dataset for each column.
This structure holds the column index and the concentration datasets for:
    NO3, NO2, SO4, pH, and EC.
"""
struct ds_m2{I, D}
    column::I
    no3::D
    #no2::D
    #so4::D
    pH::D
    ec::D
end