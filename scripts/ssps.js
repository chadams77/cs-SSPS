window.xSSPS = function(RSIZE) {

	this.RSIZE = RSIZE || 512;
	this.PCOUNT = 512;

	this.restDensity = 0.1;
	this.fieldLen = 4;
	this.gConst = 0.002;

	this.cam = {
		p: {x: 0, y: 0, z: 0},
		dir: {x: 1, y: 0, z: 0},
		up: {x: 0, y: 1, z: 0},
		nearW: 1,
		farW: 1000,
		farDist: 1000
	};
	
	this.gpu = new GPU();
	this.canvas = document.createElement('canvas');
	this.canvas.width = this.canvas.height = this.RSIZE;

	/*
		RENDER SUPPORT FUNCTIONS (GLSL)
	 */

	this.gpu.addNativeFunction('raySphere', `vec2 raySphere(vec3 r0, vec3 rd, vec3 s0, float sr) {
	    float a = dot(rd, rd);
	    vec3 s0_r0 = r0 - s0;
	    float b = 2.0 * dot(rd, s0_r0);
	    float c = dot(s0_r0, s0_r0) - (sr * sr);
	    float test = b*b - 4.0*a*c;
	    if (test < 0.0) {
	        return vec2(-1.0, 0.);
	    }
	    test = sqrt(test);
	    float X = (-b - test)/(2.0*a);
	    float Y = (-b + test)/(2.0*a);
    	return vec2(
    		X,
    		Y-X
    	);
	}`);

	this.gpu.addNativeFunction('getRay0', `vec3 getRay0(vec2 uv,  vec3 c0, vec3 cd, vec3 up, float s0, float s1, float slen) {
    	
		float _S=sin(-3.14159265359*0.5);
		float _C=cos(-3.14159265359*0.5);
		float _OC=1.0-_C;
		vec3 _AS=cd*_S;
		mat3 _P=mat3(cd.x*cd,cd.y*cd,cd.z*cd);
		mat3 _Q=mat3(_C,-_AS.z,_AS.y,_AS.z,_C,-_AS.x,-_AS.y,_AS.x,_C);
		mat3 rot = _P*_OC+_Q;

		vec3 vup = normalize(cross(cd, up));
		vec3 vleft = normalize(rot * vup);

		vec2 X = (uv - vec2(0.5, 0.5)) * s0;

		return c0 + vup * X.y + vleft * X.x;

	}`);

	this.gpu.addNativeFunction('getRayDir', `vec3 getRayDir(vec3 R0,  vec2 uv,  vec3 c0, vec3 cd, vec3 up, float s0, float s1, float slen) {
    	
		float _S=sin(-3.14159265359*0.5);
		float _C=cos(-3.14159265359*0.5);
		float _OC=1.0-_C;
		vec3 _AS=cd*_S;
		mat3 _P=mat3(cd.x*cd,cd.y*cd,cd.z*cd);
		mat3 _Q=mat3(_C,-_AS.z,_AS.y,_AS.z,_C,-_AS.x,-_AS.y,_AS.x,_C);
		mat3 rot = _P*_OC+_Q;

		vec3 vup = normalize(cross(cd, up));
		vec3 vleft = normalize(rot * vup);

		vec2 X = (uv - vec2(0.5, 0.5)) * s1;

		vec3 far = vup * X.y + vleft * X.x;

		return normalize((normalize(cd) * slen + far + c0) - R0);

	}`);

	/*
		UPDATE SUPPORT FUNCTIONS (GLSL)
	 */

	this.gpu.addNativeFunction('compGravity', `vec3 compGravity(vec3 M, vec3 J, float jMass, float f) {

		vec3 delta = J - M;
		float dlen = length(delta);
				 
		if (dlen > 0.0) {
			float dlen2 = jMass / (max(dlen*dlen, 1.) * 0.1);
			return dlen2 * (delta / dlen) * f;
		}

		return vec3(0., 0., 0.);
	}`);

	this.gpu.addNativeFunction('partDensity', `vec2 partDensity(vec3 M, vec3 J, float fLen) {

		vec3 delta = J - M;
		float len = length(delta);
				 
		if (len < fLen) {
			float t = 1. - (len / fLen);
			return vec2(t*t, t*t*t);
		}

		return vec2(0., 0.);

	}`);

	this.gpu.addNativeFunction('pressForce', `vec3 pressForce(vec3 M, vec3 J, vec3 Mv, vec3 Jv, float fLen, float dt, float spressure, float snpressure, float viscdt) {

		vec3 delta = J - M;
		float len = length(delta);		 
		if (len < fLen) {
			float t = 1. - (len / fLen);
			delta *= dt * t * (spressure + snpressure * t) / (2. * len);
			vec3 deltaV = Jv - Mv;
			deltaV *= dt * t * viscdt;
			return -(delta - deltaV);
		}

		return vec3(0., 0., 0.);

	}`);

	/*
		INIT SUPPORT FUNCTIONS
	 */	

	this.gpu.addNativeFunction('random', `float random(float sequence, float seed) {

    	return fract(sin(dot(vec2(seed, sequence), vec2(12.9898, 78.233))) * 43758.5453);
	
	}`);

	/*
		VELOCITY UPDATE KERNEL
	 */

	this.velocityKernel = this.gpu.createKernel(function(positions, velocities, attrs, dt) {

		var comp = this.thread.x % 3;
		var me = (this.thread.x - comp) / 3;
		var mx = positions[me*3],
			my = positions[me*3+1],
		    mz = positions[me*3+2];
		var mpos = [mx, my, mz];
		var mvx = velocities[me*3],
			mvy = velocities[me*3+1],
		    mvz = velocities[me*3+2];
		var mvel = [mvx, mvy, mvz];
		var mdensity = attrs[me*5+1];
		var mmass = attrs[me*5+2];
		var incomp = attrs[me*5+3];
		var viscdt = Math.pow(attrs[me*5+4], dt);
		var ddensity = 0.0,
			nddensity = 0.0;

		for (var i=0; i<this.constants.PCOUNT; i++) {
			var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
			var density = attrs[i*5+1];
			var dret = [0, 0]; dret = partDensity(mpos, opos, this.constants.fieldLen);
			ddensity += dret[0] * density;
			nddensity += dret[1] * density;
		}

		var spressure = (ddensity - this.constants.restDensity) * incomp;
		var snpressure = nddensity * incomp;

		var ret = [mvx, mvy, mvz];

		for (var i=0; i<this.constants.PCOUNT; i++) {
			if (abs(i - me) > 0.01) {
				var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
				var ovel = [velocities[i*3], velocities[i*3+1], velocities[i*3+2]];
				var mass = attrs[i*5+2];
				var jmf = (2. * mass) / (mass + mmass);
				var dret = [0., 0., 0.]; dret = pressForce(mpos, opos, mvel, ovel, this.constants.fieldLen, dt, spressure, snpressure, viscdt);
				ret[0] += dret[0]; ret[1] += dret[1]; ret[2] += dret[2];
	    	}
    	}

    	var gf = this.constants.gConst * dt;
    	for (var i=0; i<this.constants.PCOUNT; i++) {
    		var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
    		var mass = attrs[i*5+2];
    		var dret = [0., 0., 0.]; dret = compGravity(mpos, opos, mass, gf);
    		ret[0] += dret[0]; ret[1] += dret[1]; ret[2] += dret[2];
    	}

    	if (comp < 0.01) {
    		return ret[0];
    	}
    	else if (comp < 1.01) {
    		return ret[1];
    	}
    	else if (comp < 2.01) {
    		return ret[2];
    	}
	}, {
		constants: {
			PCOUNT: this.PCOUNT,
			fieldLen: this.fieldLen,
			gConst: this.gConst,
			restDensity: ((4 / 3) * Math.PI * Math.pow(this.fieldLen, 3)) * this.restDensity
		},
		output: [ this.PCOUNT*3 ],
		outputToTexture: true,
		canvas: this.canvas
	});

	/*
		POSITION UPDATE KERNEL
	 */

	this.positionKernel = this.gpu.createKernel(function(positions, velocities, dt) {
		return positions[this.thread.x] + velocities[this.thread.x] * dt;
	}, {
		output: [ this.PCOUNT*3 ],
		outputToTexture: true,
		canvas: this.canvas
	});

	/*
		RENDER KERNEL
	 */
	
	this.renderKernel = this.gpu.createKernel(function(particles, attrs, camCenter, camDir, camUp, camNearWidth, camFarWidth, camDist) {

		var uv = [this.thread.x / (this.constants.RSIZE-1), this.thread.y / (this.constants.RSIZE-1)];

		var CC = [camCenter[0], camCenter[1], camCenter[2]],
			CD = [camDir[0], camDir[1], camDir[2]],
			CU = [camUp[0], camUp[1], camUp[2]];

		var ray0 = [0,0,0]; ray0 = getRay0(uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);
		var rayDir = [0,0,0]; rayDir = getRayDir(ray0, uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);

		this.color(0., 0., 0., 1.);

		var minDist = -1.0;

		var int = 0;

		for (var i=0; i<this.constants.PCOUNT; i++) {
			var s0 = [particles[i * 3 + 0], particles[i * 3 + 1], particles[i * 3 + 2]];
			var sr = attrs[i * 5 + 0];
			if (sr > 0.0) {
				var ab = [0, 0]; ab = raySphere(ray0, rayDir, s0, sr);
				if (ab[0] >= 0.0) {
					int += (1.0 / Math.pow(1 + ab[0], 0.25))*(ab[1]/(sr*2.));
				}
			}
		}

		if (int < 0.25) {
			this.color(0., 0., int, 1.);
		}
		else if (int < 0.5) {
			this.color(0., (int-0.25)/0.25, 1., 1.);
		}
		else if (int < 1.0) {
			this.color((int-0.5)/0.5, 1., 1., 1.);
		}
		else {
			this.color(1., 1., 1., 1.);
		}

	}, {
		graphical: true,
		constants: {
			RSIZE: this.RSIZE,
			PCOUNT: this.PCOUNT,
		},
		output: [this.RSIZE, this.RSIZE],
		canvas: this.canvas,
		paramTypes: {
			camCenter: 'Array(3)',
			camDir: 'Array(3)',
			camUp: 'Array(3)',
			camNearWidth: 'Number',
			camFarWidth: 'Number',
			camDist: 'Number'
		}
	});

	/*
		INIT KERNEL
	 */

	const makeDataTexture = (array, get) => {
		const arr = [];
		for (let i=0; i<array.length; i++) {
			get(array[i], arr);
		}
		const size = arr.length;
		const kern = this.gpu.createKernel(function(parray) {
			return parray[this.thread.x];
		}, {
			outputToTexture: true,
			output: [size],
			canvas: this.canvas
		});
		return kern(arr);
	};

	const seedPositions = this.gpu.createKernel(function(){
		var t = random(Math.floor(this.thread.x / 3), this.constants.seed + 0.5);
		var r = t * this.constants.maxr;
		var a = random(Math.floor(this.thread.x / 3), this.constants.seed + 1.5);
		var comp = this.thread.x % 3;
		if (comp < 0.01) {
			return Math.cos(a) * r;
		}
		else if (comp < 1.01) {
			return Math.sin(a) * r;
		}
		else {
			return 0.0;
		}
	}, {
		constants: {
			maxr: 100,
			seed: Math.random() * 1e6
		},
		outputToTexture: true,
		output: [this.PCOUNT * 3],
		canvas: this.canvas
	});

	const seedVelocities = this.gpu.createKernel(function(){
		var t = random(this.thread.x, this.constants.seed + 0.5);
		return (t - 0.5) * this.constants.iv;
	}, {
		constants: {
			iv: 1.5,
			seed: Math.random() * 1e6
		},
		outputToTexture: true,
		output: [this.PCOUNT * 3],
		canvas: this.canvas
	});

	const seedAttrs = this.gpu.createKernel(function(){
		var comp = this.thread.x % 5;
		if (comp < 0.01) {
			return 0.25; // radius
		}
		else if (comp < 1.01) {
			return 0.5; // density
		}
		else if (comp < 2.01) {
			return 0.2; // mass
		}
		else if (comp < 3.01) {
			return 0.8; // incompress
		}
		else {
			return 0.75; // visc
		}
	}, {
		constants: {
			iv: 50,
			seed: Math.random() * 1e6
		},
		outputToTexture: true,
		output: [this.PCOUNT * 5],
		canvas: this.canvas
	});

	/*
		COPPIER
	 */
	
	const makeCoppier = (isize) => {
		return this.gpu.createKernel(function(parray) {
			return parray[this.thread.x];
		}, {
			outputToTexture: true,
			output: [isize * this.PCOUNT],
			canvas: this.canvas
		})
	};

	/*
	 * Seed Particles
	 */ 
	this.data = {
		pos: seedPositions(),
		vel: seedVelocities(),
		attr: seedAttrs()
	};

	/*
		Make coppiers
	 */
	this.posCoppier = makeCoppier(3);
	this.velCoppier = makeCoppier(3);
	this.attrCoppier = makeCoppier(5);

};


xSSPS.prototype.updateRender = function(dt) {

	// Update velocities via gravity & pressure
	this.data.vel = this.velocityKernel(
		this.data.pos,
		this.velCoppier(this.data.vel),
		this.data.attr,
		dt
	);

	// Update positions
	this.data.pos = this.positionKernel(
		this.posCoppier(this.data.pos),
		this.data.vel,
		dt
	);

	// Render
	this.renderKernel(
		this.data.pos,
		this.data.attr,
		[this.cam.p.x, this.cam.p.y, this.cam.p.z],
		[this.cam.dir.x, this.cam.dir.y, this.cam.dir.z],
		[this.cam.up.x, this.cam.up.y, this.cam.up.z],
		this.cam.nearW, this.cam.farW, this.cam.farDist
	);

};