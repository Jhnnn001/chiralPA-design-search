-- One normal-incidence Jones evaluation; lengths in µm, permittivities supplied by Python.
local function field(w)
  local k = w.k
  local norm = math.sqrt(k[1]^2+k[2]^2+k[3]^2)
  local u = w.u
  local vx = (k[2]*u[3]-k[3]*u[2])/norm
  local vy = (k[3]*u[1]-k[1]*u[3])/norm
  return {u[1]*w.cu[1]+vx*w.cv[1], u[1]*w.cu[2]+vx*w.cv[2]},
         {u[2]*w.cu[1]+vy*w.cv[1], u[2]*w.cu[2]+vy*w.cv[2]}
end

local v = {}
for token in io.read('*l'):gmatch('[^;]+') do v[#v+1] = assert(tonumber(token)) end
local S = S4.NewSimulation()
S:SetLattice({v[10],0}, {0,v[11]})
S:SetNumG(v[1])
S:AddMaterial('air', {1,0})
S:AddMaterial('TiO2', {v[13],v[14]})
S:AddMaterial('SiO2', {v[15],v[16]})
S:AddMaterial('Ag', {v[17],v[18]})
S:AddLayer('pad', 0, 'air')
S:AddLayer('inc', 0, 'air')
S:AddLayer('pat', v[8], 'air')
for i = 19, #v, 4 do
  S:SetLayerPatternRectangle('pat', 'TiO2', {v[i],v[i+1]}, 0, {v[i+2],v[i+3]})
end
S:AddLayer('sio', v[9], 'SiO2')
S:AddLayer('exit', 0, 'Ag')
S:SetFrequency(1/v[12])

local function column(s, p)
  S:SetExcitationPlanewave({0,0}, {s,0}, {p,0})
  for _, w in ipairs(S:GetWaves('inc')) do
    if w.k[3] < 0 and math.abs(w.k[1]) < 1e-6 and math.abs(w.k[2]) < 1e-6 then
      return field(w)
    end
  end
  error('No reflected specular wave')
end

local xx, yx = column(0, 1)
local xy, yy = column(1, 0)
io.write(string.format('J %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n',
  xx[1],xx[2],xy[1],xy[2],yx[1],yx[2],yy[1],yy[2]))
