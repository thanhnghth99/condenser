classdef CondenserInlet
    properties (SetAccess = immutable)
        % Declare contant input thermodynamic properties

        % R134a refrigerant side
        Refrigerant
        P_ref_in
        T_ref_in
        h_ref_in
        m_ref
        A_tube_total

        % Air side
        T_air_in
        m_air_in
        h_air
        eta_o
        A_surface_total
    end

    methods
        % Constructor
        function obj = CondenserInlet(P_ref_in, T_ref_in, m_ref, A_tube_total, T_air_in, m_air_in, h_air, eta_o, A_surface_total)
            obj.Refrigerant = 'R134a';
            obj.P_ref_in = P_ref_in;
            obj.T_ref_in = T_ref_in;
            obj.m_ref = m_ref;
            obj.A_tube_total = A_tube_total;

            obj.T_air_in = T_air_in;
            obj.m_air_in = m_air_in;
            obj.h_air = h_air;
            obj.eta_o = eta_o;
            obj.A_surface_total = A_surface_total;

            % Superheated specific enthalpy [J/kg] of inlet refrigerant to condenser
            % obj.h_ref_in = ThermoProp.get_Prop_PT('H', obj.P_ref_in, obj.T_ref_in, obj.Refrigerant);
            obj.h_ref_in = ThermoProp.get_SinglePhaseProps(obj.P_ref_in, obj.T_ref_in, obj.Refrigerant).h;
        end
    end
end