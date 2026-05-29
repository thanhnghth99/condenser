classdef CondenserModel
    properties (SetAccess = immutable)
        Inlet
        D_h
        N_c
        N_t
        A_channel
        H_cond
    end

    properties (SetAccess = protected)
        % Mass flux of refrigerant [kg/m^2.s]
        G
    end

    methods
        function obj = CondenserModel(Inlet, D_h, A_channel, N_c, N_t, H_cond)
            obj.Inlet = Inlet;
            obj.D_h = D_h;
            obj.A_channel = A_channel;
            obj.N_c = N_c;
            obj.N_t = N_t;
            obj.H_cond = H_cond;
            
            obj.G = obj.Inlet.m_ref / (A_channel * N_c * N_t);
        end

        function UA = UA(obj, w_i, h_ref_i)
            % Thermal resistance
            A_tube_total = obj.Inlet.A_tube_total;
            R_ref = 1 / (h_ref_i * A_tube_total);

            h_air = obj.Inlet.h_air;
            A_surface_total = obj.Inlet.A_surface_total;
            eta_o = obj.Inlet.eta_o;
            R_air = 1 / (h_air * A_surface_total * eta_o);

            UA = w_i / (R_ref + R_air);            
        end
    end
end