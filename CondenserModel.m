classdef CondenserModel
    properties (SetAccess = immutable)
        Inlet
        D_h
        N_c
        N_t
        W_w
        W_c
        A_channel
        W_cond
    end

    properties (SetAccess = protected)
        % Mass flux of refrigerant [kg/m^2.s]
        G
    end

    methods
        function obj = CondenserModel(Inlet, D_h, A_channel, N_c, N_t, W_w, W_c, W_cond)
            obj.Inlet = Inlet;
            obj.D_h = D_h;
            obj.A_channel = A_channel;
            obj.N_c = N_c;
            obj.N_t = N_t;
            obj.W_w = W_w;
            obj.W_c = W_c;
            obj.W_cond = W_cond;
            
            obj.G = obj.Inlet.m_ref / (A_channel * N_c * N_t);
        end

        function UA = UA(obj, w_i, h_ref_i)
            % Thermal resistance
            A_tube_total = obj.Inlet.A_tube_total;
            R_ref = 1 / (h_ref_i * A_tube_total);
            fprintf('R_ref: %.4f\n', R_ref);

            h_air = obj.Inlet.h_air;
            A_surface_total = obj.Inlet.A_surface_total;
            eta_o = obj.Inlet.eta_o;
            R_air = 1 / (h_air * A_surface_total * eta_o);
            fprintf('R_air: %.4f\n', R_air);

            UA = w_i / (R_ref + R_air);            
        end
    end
end