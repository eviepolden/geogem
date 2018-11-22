# OSBG36/EPSG27700 (British National Grid) to WGS84 (LatLong)

class NationalGrid
  attr_reader :easting, :northing

  def initialize(easting, northing)
    @easting = easting.to_i
    @northing = northing.to_i
  end

  def lat_long
    OSGB36toWGS84(easting, northing)
  end
end

def OSGB36toETRS89(x, y, h)
  #initialise constants
  e = 0.00001
  dx =
  dy =
  dz =
  x = x0-dx
  y = y0-dy
  z = z0+dz
  last_dx = dx
  last_dy = dy


end
  # def OSGB36_to_ETRS89 (x0, y0, z0 = 0.0):
  #   epsilon = 0.00001
  #   (dx, dy, dz) = _find_OSTN02_shifts_at(x0,y0)
  #   (x,  y,  z ) = (x0-dx, y0-dy, z0+dz)
  #   (last_dx, last_dy) = (dx, dy)
  #   #APPROX:
  #   while 1:
  #       (dx, dy, dz) = _find_OSTN02_shifts_at(x,y)
  #       (x, y) = (x0-dx, y0-dy)
  #       if abs(dx-last_dx)<epsilon and abs(dy-last_dy)<epsilon: break #last APPROX
  #       (last_dx, last_dy) = (dx, dy)
  #
  #   (x, y, z) = _round_to_nearest_mm(x0-dx, y0-dy, z0+dz)
  #
  #   return (x, y, z)

def ETRS89toWGS84(E,N,shape='WGS84')
  (a,b) = ellipsoid_shapes[shape]

  e2 = (a**2 - b**2) / a**2
  n = (a-b)/(a+b)

  dN = n - n0

  phi = phi0 + dN/(a * f0)

  m = compute_m(phi, b, n)
end

  # def grid_to_ll(E,N,shape='WGS84'):

	#if ( $E =~ $GR_Pattern || $E =~ $Long_GR_Pattern || $E =~ $LR_Pattern ) {
	#	($E, $N) = parse_grid($E);
	#}

	# (a,b) = ellipsoid_shapes[shape]
  #
	# e2 = (a**2.-b**2.)/a**2.
	# n = (a-b)/(a+b)
  #
	# dN = N - N0
  #
	# phi = PHI0 + dN/(a * F0)
  #
	# M = _compute_M(phi, b, n);
	# while (dN-M >= 0.001):
	#    phi = phi + (dN-M)/(a * F0)
	#    M = _compute_M(phi, b, n)
  #
	# sp2  = math.sin(phi)**2.;
	# nu   = a * F0 *			 (1. - e2 * sp2 ) ** -0.5
	# rho  = a * F0 * (1. - e2) * (1. - e2 * sp2 ) ** -1.5
	# eta2 = nu/rho - 1.
  #
	# tp = math.tan(phi)
	# tp2 = tp*tp
	# tp4 = tp2*tp2
  #
	# VII  = tp /   (2.*rho*nu)
	# VIII = tp /  (24.*rho*nu**3.) *  (5. +  3.*tp2 + eta2 - 9.*tp2*eta2)
	# IX   = tp / (720.*rho*nu**5.) * (61. + 90.*tp2 + 45.*tp4)
  #
	# sp = 1.0 / math.cos(phi)
	# tp6 = tp4*tp2
  #
	# X	= sp/nu
	# XI   = sp/(   6.*nu**3.)*(nu/rho + 2.*tp2)
	# XII  = sp/( 120.*nu**5.)*(	  5. + 28.*tp2 +   24.*tp4)
	# XIIA = sp/(5040.*nu**7.)*(	61. + 662.*tp2 + 1320.*tp4 + 720.*tp6)
  #
	# e = E - E0
  #
	# phi = phi		- VII*e**2. + VIII*e**4. -   IX*e**6.
	# lam = LAM0 + X*e -  XI*e**3. +  XII*e**5. - XIIA*e**7.
  #
	# phi = math.degrees(phi)
	# lam = math.degrees(lam)
  #
	# return (phi, lam)
	#return format_ll_ISO($phi,$lam);

def compute_m(phi, b, n)

end

def find_OSTN02_shifts(x, y)
  e_index = (x/1000)
  n_index = (y/1000)

  s0_ref = get_ostn_ref(e_index+0, n_index+0)
  s1_ref = get_ostn_ref(e_index+1, n_index+0)
  s2_ref = get_ostn_ref(e_index+0, n_index+1)
  s3_ref = get_ostn_ref(e_index+1, n_index+1)

  x0 = e_index * 1000
  y0 = e_index * 1000

  dx = x - x0
  dy = y - y0

  t = dx/1000
  u = dy/1000

  f0 = (1-t)*(1-u)
  f1 = t*(1-u)
  f2 = (1-t)*u
  f3 = t*u

  se = f0*s0_ref[0] + f1*s1_ref[0] + f2*s2_ref[0] + f3*s3_ref[0]
  sn = f0*s0_ref[1] + f1*s1_ref[1] + f2*s2_ref[1] + f3*s3_ref[1]
  sg = f0*s0_ref[2] + f1*s1_ref[2] + f2*s2_ref[2] + f3*s3_ref[2]

  return (se, sn, sg)
end

def get_ostn_ref(x,y)
  key = "#{y}03x#{x}03x"
  # key_data = (easting / 1000)
  # key_data = (data[0]/1000.0 + MIN_X_SHIFT,data[1]/1000.0 +MIN_Y_SHIFT,data[2]/1000.0 + MIN_Z_SHIFT)
  # ostn_shift_for[key] = key_data
  # return key_data
end

  # def _get_ostn_ref(x,y):
  #
  #     key = "%03x%03x" % (y, x)
  #     if key in ostn_shift_for:
  #         return ostn_shift_for[key]
  #
  #     if key in ostn_data:
  #         data = ostn_data[key]
  #         data2 = (data[0]/1000.0 + MIN_X_SHIFT,data[1]/1000.0 +MIN_Y_SHIFT,data[2]/1000.0 + MIN_Z_SHIFT)
  #         ostn_shift_for[key] = data2
  #         return data2


def OSGB36toWGS84(easting, northing)
  #Initialising ellipsoidal and projection constants
  airy_semimajor = 6377563.396 #Airy 1830 semi-major axis used for OSGB36 (m)
  airy_semiminor = 6356256.909 #Airy 1830 semi-minor axis used for OSGB36 (m)
  scale_factor = 0.9996012717 #scale factor on the central meridian
  pi = Math::PI
  lat_origin = (49*pi)/180 #latitude of true origin (radians)
  lon_origin = (-2*pi)/180 #longtitude of true origin (radians)
  northing_origin = -100000 #northing of true origin (m)
  easting_origin = 400000 #easting of true origin (m)

  e_squared = 1 - ((airy_semiminor * airy_semiminor) / (airy_semimajor * airy_semimajor)) #eccentricity squared
  n = (airy_semimajor - airy_semiminor) / (airy_semimajor + airy_semiminor)

  #iteration
  lat = lat_origin
  m = 0

  lat = (n-northing_origin-m)/(airy_semimajor*scale_factor) + lat
  m_1 = (1 + n + (5./4)*n**2 + (5./4)*n**3) * (lat-lat_origin)
  m_2 = (3*n + 3*n**2 + (21./8)*n**3) * Math.sin(lat-lat_origin) * Math.cos(lat+lat_origin)
  m_3 = ((15./8)*n**2 + (15./8)*n**3) * Math.sin(2*(lat-lat_origin)) * Math.cos(2*(lat+lat_origin))
  m_4 = (35./24)*n**3 * Math.sin(3*(lat-lat_origin)) * Math.cos(3*(lat+lat_origin))

  # meridional arc
  m = airy_semiminor * scale_factor * (m_1 - m_2 + m_3 - m_4)

  #transverse radius of curvature
  nu = (airy_semimajor*scale_factor)/Math.sqrt(1-e_squared*Math.sin(lat)**2)

  #meridional radius of curvature
  rho = airy_semimajor*scale_factor*(1-e_squared)*(1-e_squared*Math.sin(lat)**2)**(-1.5)
  eta_squared = nu/rho-1

  sec_lat = 1/Math.cos(lat)
  vii = Math.tan(lat)/(2*rho*nu)
  viii = Math.tan(lat)/(24*rho*nu**3)*(5+3*Math.tan(lat)**2+eta_squared-9*Math.tan(lat)**2*eta_squared)
  ix = Math.tan(lat)/(720*rho*nu**5)*(61+90*Math.tan(lat)**2+45*Math.tan(lat)**4)
  x = sec_lat/nu
  xi = sec_lat/(6*nu**3)*(nu/rho+2+Math.tan(lat)**2)
  xii = sec_lat/(120*nu**5)*(5+28*Math.tan(lat)**2+24*Math.tan(lat)**4)
  xiia = sec_lat/(5040*nu**7)*(61+662*Math.tan(lat)**2+1320*Math.tan(lat)**4+720*Math.tan(lat)**6)
  dE = easting - easting_origin

  #Airy 1830 (denoted by _1)
  lat_1 = lat - vii*dE**2 + viii*dE**4 - ix*dE**6
  lon_1 = lon_origin + x*dE - xi*dE**3 + xii*dE**5 - xiia*dE**7

  #conversion to GRS80 ellipsoid
  #conversion to cartesian coordinates
  h = 0 #third cartesian coord
  x_1 = (nu/scale_factor + h)*Math.cos(lat_1)*Math.cos(lon_1)
  y_1 = (nu/scale_factor + h)*Math.cos(lat_1)*Math.sin(lon_1)
  z_1 = ((1-e_squared)*nu/scale_factor + h)*Math.sin(lat_1)

  #Helmert transformation to go from Airy 1830 -> GRS80
  s = -20.4894*10**-6 #scale factor
  tx = 446.448
  ty = -125.157
  tz = 542.060
  rxs = 0.1502
  rys = 0.2470
  rzs = 0.8421

  #convert seconds to radians
  def sec_to_rad(x)
    x*Math::PI/(180*3600)
  end

  rx = sec_to_rad(rxs)
  ry = sec_to_rad(rys)
  rz = sec_to_rad(rzs)
  x_2 = tx + (1+s)*x_1 + (-rz)*y_1 + (ry)*z_1
  y_2 = ty + (rz)*x_1 + (1+s)*y_1 + (-rx)*z_1
  z_2 = tz + (-ry)*x_1 + (rx)*y_1 + (1+s)*z_1

  #back to polar coordinates
  a_2 = 6378137.000
  b_2 = 6356752.3141
  e_squared_2 = 1 - (b_2*b_2)/(a_2*a_2) # The eccentricity(squared) of the GRS80 ellipsoid
  p_1 = Math.sqrt(x_2**2 + y_2**2)

  #obtain latitude
  lat = Math.atan2(z_2, (p_1*(1-e_squared_2)))
  latold = 2*pi

  nu_2 = a_2/Math.sqrt(1-e_squared_2*Math.sin(latold)**2)
  lat = Math.atan2(z_2+e_squared_2*nu_2*Math.sin(latold), p_1)

  #obtain longitude and height
  lon = Math.atan2(y_2, x_2)
  h = p_1/Math.cos(lat) - nu_2

  #convert to degrees
  lat = lat*180/pi
  lon = lon*180/pi

  return [lat, lon]
end
