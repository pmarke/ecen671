classdef rlsDynamics < handle
    %----------------------------
    properties
        state
        R_inv
        k
        A
        c
        y_history
        Ts
        x     % [-1/ao; -a1/ao; bo/ao]
        f     % [yddot, ydot, u]
        
        
        y     % y[i]
        ydl   % y[i-1]
        yddl  % y[i-2]
        ydddl
        ydot_dl
        udl
        uddl
        
    end
    %----------------------------
    methods
        %---constructor-------------------------
        function self = rlsDynamics(Ts)
            % Initial state conditions
            self.R_inv = eye(3)*1e-10;
            self.A = [];
            self.y_history = [];
            self.c = [0,0,0]'; 
            self.Ts = Ts; 
            self.state = [0;0;0];
            self.x = [0;0;0];
            self.f = [0,0,0];
            self.y = 0;
            self.ydl = 0;
            self.yddl = 0;
            self.ydddl = 0;
            self.ydot_dl = 0;
            self.udl = 0;
            self.uddl = 0;
          
        end
        %----------------------------
        function self = propagateDynamics(self, u,f)
            %
            % Use difference equations to construct f
            % 
            
            ydot =  (self.y - self.ydl)/self.Ts;
            yddot = (self.y -2*self.ydl + self.yddl)/self.Ts^2;
            self.f = [yddot, self.ydot_dl, self.uddl];

            
            % Update terms
            self.ydddl = self.yddl;
            self.yddl = self.ydl;
            self.ydl = self.y;            
            self.ydot_dl = ydot;
            self.uddl = self.udl;
            self.udl = u;
            
            
        end
        %----------------------------
        function self = update(self, y, u,f)
            
            self.y = y;
            self.propagateDynamics(u,f);  % Construct f
            
            % Not enough data points to do RLS
            if(size(self.A,1) < 5)
                self.y_history = [self.ydddl;self.y_history];
                self.A = [self.f;self.A];
                self.R_inv = inv(self.A'*self.A + 1e-10*eye(3));
                self.c = self.A'*self.y_history;
                self.x = self.R_inv*self.c;
            
            % You have enough data points so that R_inv is non_zero
            else
                self.k = self.R_inv*self.f'/(1+self.f*self.R_inv*self.f');
                self.R_inv = self.R_inv - self.k*self.f*self.R_inv;
%                 R_inv = self.R_inv
                self.c = self.c + self.f'*self.ydddl;
            
                self.x = self.R_inv*self.c;
            end
            
            

            
            
        end
        %----------------------------
        function c = states(self)
            
            a0 = -1/self.x(1);
            a1 = self.x(2)/(self.x(1)+1e-10);
            b0 =  self.x(3)/(-self.x(1)+1e-10);
            
            c = [a1,a0,b0];
            
        end
    end
end