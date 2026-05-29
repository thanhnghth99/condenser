classdef ThermoProp
    % Utility class specifically to call the CoolProp library for thermodynamics properties

    methods (Static)
        % Calculate saturated temperature [K]
        function T_sat = get_T_sat(P, x, fluid)
            T_sat = py.CoolProp.CoolProp.PropsSI('T', 'P', P, 'Q', x, fluid);            
        end

        % Calculate saturated specific enthalpy [J/kg]
        function h_sat = get_h_sat(P, x, fluid)
            h_sat = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', x, fluid);            
        end

        % Calculate single phase properties
        function props = get_SinglePhaseProps(P, T, fluid)
            
            % Constant pressure specific heat capacity[J/kg.K]
            props.cp = py.CoolProp.CoolProp.PropsSI('C', 'P', P, 'T', T, fluid);

            % Viscosity [Pa.s]
            props.mu = py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'T', T, fluid);

            % Thermal conductivity [W/m.K]
            props.k = py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'T', T, fluid);

            % Mass density [kg/m^3]
            props.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'T', T, fluid);

            % Specific enthalpy [J/kg]
            props.h = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'T', T, fluid);

            % Prandtl number
            props.Pr = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P, 'T', T, fluid);
        end

        % Extract the desired properties
        function value = get_Prop_PT(desired_prop, P, T, fluid)
            value = py.CoolProp.CoolProp.PropsSI(desired_prop, 'P', P, 'T', T, fluid);
        end
    end
end