function linkages(scene)
if nargin < 1
    scene = 4;
end

% Set up the scene here.
% Note that links don't have to be placed exactly. The first call to
% solveLinkage() will automatically snap the links so that all the
% constraints are satisfied.
linkages = [];
pins = [];
particles = [];
switch scene
    case 0
        % Crank-rocker with a fixed range
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2]; % link between link 1 and 2
        pins(1).pts = [0,0;0,0]'; % transposed! link 1 (0, 0) connects with link 2 (0, 0)
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [1+rand(1),0;2,0]'; % pin location on link3 is randomized
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 4; % which link?
        particles(1).pt = [0.5;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
        particles(2).link = 4;
        particles(2).pt = [2.5;-0.1];
        particles(2).ptsWorld = zeros(2,0);
        
        isCounterClockwise = 1;
        
        % Angle range 
        angStart = -pi/2;
        angEnd = pi/2;
        % driver angular velocity
        angVel = 2*pi;
    case 1
        % Drag-link mech
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;0.4,-0.1;0.4,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/2;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;1.2,-0.1;1.2,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;1.4,-0.1;1.4,0.1;0,0.1]';
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2]; % link between link 1 and 2
        pins(1).pts = [0,0;0,0]'; % transposed! link 1 (0, 0) connects with link 2 (0, 0)
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [0.4,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [1.2,0;1.4,0]'; 
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 4; % which link?
        particles(1).pt = [1.4;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
        particles(2).link = 4;
        particles(2).pt = [0;0.1];
        particles(2).ptsWorld = zeros(2,0);
        
        isCounterClockwise = 1;
        
        % Angle range 
        angStart = intmin;
        angEnd = intmax;
        % driver angular velocity
        angVel = 2*pi;
    case 2
        % Double-rocker
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/4;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/3;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2]; % link between link 1 and 2
        pins(1).pts = [0,0;0,0]'; % transposed! link 1 (0, 0) connects with link 2 (0, 0)
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [1,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [3,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [2.5,0;1,0]'; 
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 3; % which link?
        particles(1).pt = [2.5;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
        particles(2).link = 2;
        particles(2).pt = [3;0.1];
        particles(2).ptsWorld = zeros(2,0);
        
        isCounterClockwise = 1;
        
        % Angle range 
        angStart = pi/4;
        angEnd = pi/2;
        % driver angular velocity
        angVel = 2*pi;
    case 3
        % Hoekens
        % Bottom link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;2,-0.1;2,0.1;0,0.1]'; % display vertices
        % Left link
        links(2).angle = pi/2;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;1,-0.1;1,0.1;0,0.1]';
        % Right link
        links(3).angle = pi/3;
        links(3).pos = [1;0];
        links(3).verts = [0,-0.1;2.5,-0.1;2.5,0.1;0,0.1]';
        % Top link
        links(4).angle = 0;
        links(4).pos = [-1;1];
        links(4).verts = [0,-0.1;5,-0.1;5,0.1;0,0.1]';
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2]; % link between link 1 and 2
        pins(1).pts = [0,0;0,0]'; % transposed! link 1 (0, 0) connects with link 2 (0, 0)
        % Bottom-right
        pins(2).links = [1,3];
        pins(2).pts = [2,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [1,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [2.5,0;2.5,0]'; 
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 2; % which link?
        particles(1).pt = [1;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
        particles(2).link = 4;
        particles(2).pt = [5;0];
        particles(2).ptsWorld = zeros(2,0);
        
        isCounterClockwise = 1;
        
        % Angle range 
        angStart = intmin;
        angEnd = intmax;
        % driver angular velocity
        angVel = 2*pi;
    case 4
        % Peaucellier-Lipkin
        % Bottom ground link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]'; % display vertices
        % drive link
        links(2).angle = 0;
        links(2).pos = [2;0];
        links(2).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        

        % Lower left link
        links(3).angle = -pi/4;
        links(3).pos = [5;0];
        links(3).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Lower right link
        links(4).angle = pi/4;
        links(4).pos = [5+3/sqrt(2);-3/sqrt(2)];
        links(4).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Upper left link
        links(5).angle = pi/4;
        links(5).pos = [5;0];
        links(5).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
        % Upper right link
        links(6).angle = 7/4*pi;
        links(6).pos = [5+3/sqrt(2);3/sqrt(2)];
        links(6).verts = [0,-0.1;3,-0.1;3,0.1;0,0.1]';
       
        % Upper Support link
        links(7).angle = pi/6;
        links(7).pos = [-1;0];
        sup_length = sqrt((3+3+3/sqrt(2))^2 + (3/sqrt(2))^2);
        links(7).verts = [0,-0.1;sup_length,-0.1;sup_length,0.1;0,0.1]';
        % Lower Support link
         links(8).angle = -pi/6
         links(8).pos = [-1;0];
         links(8).verts = [0,-0.1;sup_length,-0.1;sup_length,0.1;0,0.1]';
        
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % 1-2
        pins(1).links = [1,2]; % link between link 1 and 2
        pins(1).pts = [3,0;0,0]'; % transposed! link 1 (0, 0) connects with link 2 (0, 0)
        % 1-7
        pins(2).links = [1,7];
        pins(2).pts = [0,0;0,0]';
        % 1-8
        pins(3).links = [1,8];
        pins(3).pts = [0,0;0,0]';
         % 2-3
        pins(4).links = [2,3];
        pins(4).pts = [3,0;0,0]';
         % 2-5
        pins(5).links = [2,5];
        pins(5).pts = [3,0;0,0]';
         % 3-4
        pins(6).links = [3,4];
        pins(6).pts = [3,0;0,0]';
         % 3-8
        pins(7).links = [3,8];
        pins(7).pts = [3,0;sup_length,0]';
         % 5-6
        pins(8).links = [5,6];
        pins(8).pts = [3,0;0,0]';
         % 6-4
        pins(9).links = [6,4];
        pins(9).pts = [3,0;3,0]';
         % 7-5
        pins(10).links = [7,5];
        pins(10).pts = [sup_length,0;3,0]';
         % 8-3
        pins(11).links = [8,3];
        pins(11).pts = [sup_length,0;3,0]';
         % 4-8
        pins(12).links = [4,8];
        pins(12).pts = [0,0;sup_length,0]';
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 6; % which link?
        particles(1).pt = [3;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        
        particles(2).link = 2;
        particles(2).pt = [3;0];
        particles(2).ptsWorld = zeros(2,0);
        
        isCounterClockwise = 1;
        
        % Angle range 
        angStart = -pi/4;
        angEnd = pi/4;
        % driver angular velocity
        angVel = 2*pi;
    case 5
        % Klann
        % Ground link
        links(1).angle = 0; % rotation from the positive x-axis
        links(1).pos = [27.690;31.092]; % position of the center of rotation
        links(1).verts = [0,0;0,-11.387;13.705,-4.451]'; % display vertices
        % Driver link
        links(2).angle = pi;
        links(2).pos = [41.395;26.641];
        links(2).verts = [0,-0.1;5.500,-0.1;5.500,0.1;0,0.1]';
        % Middle triangle
        links(3).angle = 0;
        links(3).pos = [36.033;25.415];
        links(3).verts = [0,0;-25.154,-3.497;-14.532,-3.326]';
        % Lower link
        links(4).angle = pi;
        links(4).pos = [27.690;19.705];
        links(4).verts = [0,-0.1;6.63,-0.1;6.63,0.1;0,0.1]';
        % Upper link
        links(5).angle = pi;
        links(5).pos = [27.690;31.092];
        links(5).verts = [0,-0.1;8.45,-0.1;8.45,0.1;0,0.1]';
        % Left triangle
        links(6).angle = 0;
        links(6).pos = [19.432;32.923];
        links(6).verts = [0,0;-8.55,-11.005;-19.432,-32.923]';
        
        % Which link is grounded? - 1
        grounded = 1;
        % Which link is the driver? - 2
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % 1-2
        pins(1).links = [1,2]; 
        pins(1).pts = [13.705,-4.451;0,0]';
        % 2-3
        pins(2).links = [2,3]; 
        pins(2).pts = [5.5,0;0,0]';
        % 4-3
        pins(3).links = [4,3]; 
        pins(3).pts = [6.63,0;-14.532,-3.326]';
        % 1-5
        pins(4).links = [1,5]; 
        pins(4).pts = [0,0;0,0]';
        % 5-6
        pins(5).links = [5,6]; 
        pins(5).pts = [8.45,0;0,0]';
        % 3-6
        pins(6).links = [3,6]; 
        pins(6).pts = [-25.154,-3.497;-8.55,-11.005]';
        % 1-4
        pins(7).links = [1,4]; 
        pins(7).pts = [0,-11.387;0,0]';
        
        % List of tracer particles for display. Used to draw the track of
        % the links
        particles(1).link = 2; % which link?
        particles(1).pt = [5.5;0]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
       
        particles(2).link = 6; % which link?
        particles(2).pt = [-19.432;-32.923]; % tracer particle point in local coords
        particles(2).ptsWorld = zeros(2,0); % transformed points, initially empty
       
        isCounterClockwise = -1;
        
        % Angle range 
        angStart = intmin;
        angEnd = intmax;
        % driver angular velocity
        angVel = 2*pi;
        
    case 10
        % Extra credit for scissor mechanism
        % Ground link
       
end

% Initialize
for i = 1 : length(links)
    links(i).grounded = (i == grounded);
    links(i).driver = (i == driver);
    % These target quantities are only used for grounded and driver links
    links(i).angleTarget = links(i).angle;
    links(i).posTarget = links(i).pos;
end
drawScene(links,pins,particles,true);

% lsqnonlin options
if verLessThan('matlab','8.1')
    opt = optimset(...
        'Jacobian','on',... % off for scissors
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
else
    opt = optimoptions('lsqnonlin',...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
end

% Simulation loop
t = 0; % current time
T = 1; % final time
dt = 0.01; % time step
while t < T
    % Procedurally set the driver angle.
    % Right now, the target angle is being linearly increased, but you may
    % want to do something else.
    links(driver).angleTarget = links(driver).angleTarget + dt*angVel*isCounterClockwise;
    if links(driver).angleTarget <= angStart & isCounterClockwise == -1
        isCounterClockwise = 1;
    elseif links(driver).angleTarget >= angEnd & isCounterClockwise == 1
        isCounterClockwise = -1;
    end
    % Solve for linkage orientations and positions
    [links,feasible] = solveLinkage(links,pins,opt); % for extra credit?
    % Update particle positions
    particles = updateParticles(links,particles);
    % Draw scene
    drawScene(links,pins,particles);
    % Quit if over-constrained
    if ~feasible
        break;
    end
    t = t + dt;
end


% 
% Simulation loop
% t = 0; % current time
% T = 1; % final time
% dt = 0.01; % time step
% angVel = 2*pi; % driver angular velocity
% while t < T
%     Procedurally set the driver angle.
%     Right now, the target angle is being linearly increased, but you may
%     want to do something else.
%     links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
%     if links(driver).angleTarget == pi / 2
%         
%     Solve for linkage orientations and positions
%     [links,feasible] = solveLinkage(links,pins,opt); % for extra credit!!!!!!!!!!!!!!!!
%     Update particle positions
%     particles = updateParticles(links,particles);
%     Draw scene
%     drawScene(links,pins,particles);
%     Quit if over-constrained
%     if ~feasible
%         break;
%     end
%     t = t + dt;
% end
    

end

%%
function [R,dR] = rotationMatrix(angle) % return R and dR
c = cos(angle);
s = sin(angle);
% Rotation matrix
R = zeros(2);
R(1,1) = c;
R(1,2) = -s;
R(2,1) = s;
R(2,2) = c;
if nargout >= 2
    % Rotation matrix derivative
    dR = zeros(2);
    dR(1,1) = -s;
    dR(1,2) = -c;
    dR(2,1) = c;
    dR(2,2) = -s;
end
end

%%
function [links,feasible] = solveLinkage(links,pins,opt)
nlinks = length(links);
% Extract the current angles and positions into a vector
angPos0 = zeros(3*nlinks,1);
for i = 1 : nlinks
    link = links(i);
    ii = (i-1)*3+(1:3);
    angPos0(ii(1)) = link.angle;
    angPos0(ii(2:3)) = link.pos;
end
% Limits
lb = -inf(size(angPos0));
ub =  inf(size(angPos0));
% Solve for angles and positions
[angPos,r2] = lsqnonlin(@(angPos)objFun(angPos,links,pins),angPos0,lb,ub,opt);
% If the mechanism is feasible, then the residual should be zero
feasible = true;
if r2 > 1e-6 % it cannot be larger than 1. Maybe a slightly larger value is enough
    fprintf('Mechanism is over constrained!\n');
    feasible = false;
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end
end

%% needed for scissors.
function [f,J] = objFun(angPos,links,pins)
nlinks = length(links);
npins = length(pins);
% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins; % 3 for ground, 3 for driver, 2*npins for pins
f = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
    link = links(i);
    if link.grounded || link.driver
        % Grounded and driver links have their angles and positions
        % prescribed.
        f(k+1,    1) = link.angle - link.angleTarget;
        f(k+(2:3),1) = link.pos - link.posTarget;
        % The Jacobian of this constraint is the identity matrix
        J(k+(1:3),k+(1:3)) = eye(3);
        k = k + 3;
    end
end
% Pin constraints
for i = 1 : npins
    pin = pins(i);
    rows = k+(1:2); % row index of this pin constraint
    indLinkA = pin.links(1); % array index of link A
    indLinkB = pin.links(2); % array index of link B
    linkA = links(indLinkA);
    linkB = links(indLinkB);
    [Ra,dRa] = rotationMatrix(linkA.angle);
    [Rb,dRb] = rotationMatrix(linkB.angle);
    % Local positions
    ra = pin.pts(:,1);
    rb = pin.pts(:,2);
    % World positions
    xa = Ra * ra + linkA.pos;
    xb = Rb * rb + linkB.pos;
    p = xa(1:2) - xb(1:2);
    f(rows,1) = p;
    % Column indices for the angles and positions of links A and B
    colAngA = (indLinkA-1)*3 + 1;
    colPosA = (indLinkA-1)*3 + (2:3);
    colAngB = (indLinkB-1)*3 + 1;
    colPosB = (indLinkB-1)*3 + (2:3);
    % The Jacobian of this constraint is the partial derivative of f wrt
    % the angles and positions of the two links.
    J(rows,colAngA) = dRa * ra;
    J(rows,colPosA) = eye(2);
    J(rows,colAngB) = -dRb * rb;
    J(rows,colPosB) = -eye(2);
    k = k + 2;
end
end

%%
function particles = updateParticles(links,particles)
% Transform particle position from local to world
for i = 1 : length(particles)
    particle = particles(i);
    link = links(particle.link);
    R = rotationMatrix(link.angle);
    x = R * particle.pt + link.pos;
    % Append world position to the array (grows indefinitely)
    particles(i).ptsWorld(:,end+1) = x;
end
end

%%
function drawScene(links,pins,particles,initialize)
if nargin < 4
    initialize = false;
end
if initialize
    clf;
else
    cla;
end
hold on;
grid on;
% Draw links
for i = 1 : length(links)
    link = links(i);
    R = rotationMatrix(link.angle);
    % Draw frame
    p = link.pos; % frame origin
    s = 0.1; % frame display size
    px = p + s*R(:,1); % frame x-axis
    py = p + s*R(:,2); % frame y-axis
    plot([p(1),px(1)],[p(2),px(2)],'r','LineWidth',3);
    plot([p(1),py(1)],[p(2),py(2)],'g','LineWidth',3);
    % Draw link geometry
    if link.grounded
        color = [1 0 0];
    elseif link.driver
        color = [0 1 0];
    else
        color = [0 0 1];
    end
    E = [R,link.pos;0,0,1]; % transformation matrix
    vertsLocal = [link.verts;ones(1,size(link.verts,2))];
    vertsWorld = E * vertsLocal;
    plot(vertsWorld(1,[1:end,1]),vertsWorld(2,[1:end,1]),'Color',color);
end
% Draw pins
for i = 1 : length(pins)
    pin = pins(i);
    linkA = links(pin.links(1));
    linkB = links(pin.links(2));
    Ra = rotationMatrix(linkA.angle);
    Rb = rotationMatrix(linkB.angle);
    xa = Ra * pin.pts(:,1) + linkA.pos;
    xb = Rb * pin.pts(:,2) + linkB.pos;
    plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','c');
    plot(xb(1),xb(2),'mx','MarkerSize',10,'LineWidth',2);
end
% Draw particles
for i = 1 : length(particles)
    particle = particles(i);
    if ~isempty(particle.ptsWorld)
        plot(particle.ptsWorld(1,:),particle.ptsWorld(2,:),'k');
        plot(particle.ptsWorld(1,end),particle.ptsWorld(2,end),'ko');
    end
end
axis equal;
drawnow;
end
