within gantry_system;

model trolley_pendulum
  type Length = Real(unit="m");
  type Velocity = Real(unit="m/s");
  type AngularDisp = Real(unit="rad");
  type AngularVelo = Real(unit="rad/s");
  type Mass = Real(unit="kg", min=0);
  type Damping = Real;
  type Acceleration = Real(unit="m/s2");
  type ControlSignal = Integer;

  parameter Mass m=200 "Mass of pendulum bob/container";
  parameter Mass M=10 "Mass of trolley/cart";
  parameter Length r=1 "Length of the rope connecting the pendulum bob to the trolley";
  parameter Damping d_p=0.5 "Damping factor swinging of pendulum";
  parameter Damping d_c=2 "Damping factor for motion of cart";
  parameter Acceleration g=9.8 "Constant for gravitational acceleration on the surface of the Earth (not Newton's gravitational constant!!)";
  
  // Variables
  Length x=0 "Displacement of the trolley/cart";
  Velocity v=0 "Velocity of the trolley/cart";
  AngularDisp tetha=0 "Angular displacement of the pendulum, w.r.t the trolley";
  AngularVelo omega=0 "Angular velocity of the pendulum";
  ControlSignal u=0 "Control signal to move the trolley and pendulum";
equation
  der(x) = v;
  der(tetha) = omega;
  der(v) = ( r*(d_c*v - m*(g*sin(tetha)*cos(tetha) + r*sin(tetha)*omega^2) - u) - (d_p*cos(tetha)*omega) ) / ( -r*(M+m*sin(tetha)^2) );
  der(omega) = ( (d_p*omega*(m+M)) + (m^2*r^2*sin(tetha)*cos(tetha)*omega^2) + m*r*( (g*sin(tetha)*(m+M) + ( cos(tetha)*(u-d_c*v)))) ) / ( (m*r^2)*(-M-(m*sin(tetha)^2)) );

end trolley_pendulum;