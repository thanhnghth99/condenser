classdef ThermoProp
    % Utility class specifically to call the CoolProp library for thermodynamics properties

    methods (Static)
        % Calculate saturated temperature [K]
        function T_sat = get_T_sat(P, x, fluid)
            T_sat = py.CoolProp.CoolProp.PropsSI('T', 'P', P, 'Q', x, fluid);
        end

        % Calculate saturated vapor thermodynamic properties
        function vapor_props = get_SatVaporProps(P, fluid)
            
            % Viscosity [Pa.s]
            vapor_props.mu_v = py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'Q', 1, fluid);

            % Thermal conductivity [W/m.K]
            vapor_props.k_v = py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'Q', 1, fluid);

            % Mass density [kg/m^3]
            vapor_props.rho_v = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 1, fluid);

            % Specific enthalpy [J/kg]
            vapor_props.h_v = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 1, fluid);

            % Prandtl number
            vapor_props.Pr_v = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P, 'Q', 1, fluid);
        end

        % Calculate saturated liquid thermodynamic properties
        function liquid_props = get_SatLiquidProps(P, fluid)
            
            % Viscosity [Pa.s]
            liquid_props.mu_l = py.CoolProp.CoolProp.PropsSI('V', 'P', P, 'Q', 0, fluid);

            % Thermal conductivity [W/m.K]
            liquid_props.k_l = py.CoolProp.CoolProp.PropsSI('L', 'P', P, 'Q', 0, fluid);

            % Mass density [kg/m^3]
            liquid_props.rho_l = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 0, fluid);

            % Specific enthalpy [J/kg]
            liquid_props.h_l = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 0, fluid);

            % Prandtl number
            liquid_props.Pr_l = py.CoolProp.CoolProp.PropsSI('PRANDTL', 'P', P, 'Q', 0, fluid);
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