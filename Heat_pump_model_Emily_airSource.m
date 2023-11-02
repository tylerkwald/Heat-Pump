function [COP,eta2,Sdot_gen_cond,Sdot_gen_valve,Sdot_gen_comp,Sdot_gen_loop,Sdot_gen_whole] = Heat_pump_model_Emily(fluid,m_dot_r,eta,T_pinchC,T_pinchE,T_g,T_i,T_w,T_r,T_h)
%Emily's non-idealized heat pump COP function, updated 2/20/23

%define evaporation/condensation temperatures and pressures
T_L=T_r-T_pinchE;%bottom temperature of HP (temp in evap)
T_H=T_h+T_pinchC;%top saturation temperature of HP (excluding superheat) (temp in cond)
P_L=refpropm('P','T',T_L,'Q',1,fluid);%pressure (in kPa) of evaporator
P_H=refpropm('P','T',T_H,'Q',1,fluid);%pressure (in kPa) of condenser

%Use components to sequentially find state enthalpies for COP
[h_in_comp,h_out_comp,T_SH,Sdot_gen_comp] = realcompressor(fluid,P_L,P_H,eta,m_dot_r);
[Sdot_gen_cond,h_out_cond,m_dot_a,c_pa] = condenser_airheater(fluid,P_H,T_SH,m_dot_r,T_i,T_h);
%calculate COP from enthalpies using Equation 2.5 in Sarbu and Sebarchievici
COP=(h_out_comp-h_out_cond)/(h_out_comp-h_in_comp);

eta2=COP*(T_i-T_g)/T_i;

%calculate entropy generation (for a and c see functions at end)
%b. valve
s_vi=refpropm('S','T',T_H,'Q',0,fluid);%specific entropy in (J/kg-K)
s_vo=refpropm('S','P',P_L,'H',h_out_cond,fluid);%specific entropy out (J/kg-K)
Sdot_gen_valve=m_dot_r*(s_vo-s_vi);

%d. ground loop
h_in=h_out_cond;%enthalpy of refrig stream at evap (and valve) inlet (in J/kg)
h_out=refpropm('H','T',T_L,'Q',1,fluid);%enthalpy of refrig stream at evap inlet (in J/kg)
h_w_in=refpropm('H','T',T_w,'P',100,'air.mix');%enthalpy of air at evap inlet
h_w_out=refpropm('H','T',T_r,'P',100,'air.mix');%enthalpy of air at evap outlet
Q_dotg = m_dot_r*(h_out-h_in);%heat transfer rate into refrigerant and into ground loop from ground
m_dot_w=Q_dotg/(h_w_in-h_w_out);%air mass flow rate
s_w_intoground=refpropm('S','T',T_r,'P',100,'air.mix');
s_w_outofground=refpropm('S','T',T_w,'P',100,'air.mix');
Sdot_gen_loop=m_dot_w*(s_w_outofground-s_w_intoground)-Q_dotg/T_g;

%e. whole thing
%Sdot_gen_whole=m_dot_a*c_pa*log(T_h/T_i)-Q_dotg/T_g;
Sdot_gen_whole=m_dot_r*(h_out_comp-h_in_comp)/T_g+m_dot_a*c_pa*(T_i-T_h)/T_g+m_dot_a*c_pa*log(T_h/T_i);%how to get ground heat transfer rate if you skipped part d

end

function [h_in,h_out,T_out,Sdot_gen_comp] = realcompressor(fluid,P_1,P_2,eta,m_dot_r)
%define inlet state
%fluid enters saturated
h_in=refpropm('H','P',P_1,'Q',1,fluid);
s_in=refpropm('S','P',P_1,'Q',1,fluid);

%hypothetical outlet state of reversible turbine
h_out_rev=refpropm('H','P',P_2,'S',s_in,fluid);

%real outlet
h_out=(h_out_rev-h_in)/eta+h_in;
T_out=refpropm('T','P',P_2,'H',h_out,fluid);
s_out=refpropm('S','P',P_2,'H',h_out,fluid);
Sdot_gen_comp=m_dot_r*(s_out-s_in);
end

function [Sdot_gen_cond,h_out,m_dot_a,c_pa] = condenser_airheater(fluid,P_r,T_rin,m_dot_r,T_a_in,T_a_out)
%Calculate properties
%inlet
h_in=refpropm('H','T',T_rin,'P',P_r,fluid);%enthalpy of refrig stream at inlet (in J/kg)
%outlet
h_out=refpropm('H','P',P_r,'Q',0,fluid);%enthalpy of refrig stream at outlet (in J/kg)
%general air
c_pa=1000;%specific heat of air at average temperature (in J/kg-K)

%calculate remaining outputs
Q_dot = m_dot_r*(h_out-h_in);%heat transfer rate *into* refrigerant (in W; will be <0)
m_dot_a=Q_dot/(c_pa*(T_a_in-T_a_out));%air mass flow rate (in kg/s)

%entropy
s_r_in=refpropm('S','T',T_rin,'P',P_r,fluid);%entropy of refrig stream at inlet (in J/kg-K)
%outlet
s_r_out=refpropm('S','P',P_r,'Q',0,fluid);%entrolpy of refrig stream at outlet (in J/kg-K)
%general air
Sdot_gen_cond=m_dot_a*c_pa*log(T_a_out/T_a_in)+m_dot_r*(s_r_out-s_r_in);
end