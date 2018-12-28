window.xSSPS = function(scene, camera) {

	this.list = [];
	this.pgeom2 = new THREE.IcosahedronGeometry(0.5, 0);
	this.pgeom = new THREE.PlaneBufferGeometry(3, 3, 1, 1);
	this.hSize = 0.75;
	this.mhSize = 2.5;
	this.nGravity = 15;
	this.gConst = 0.005;
	this.lFactor = 20.0;
	this.minMassGravity = 10;
	this.scene = scene;
	this.camera = camera;

	this.seedTypes();

};

const sphereDist = function(ax, ay, az, bx, by, bz, ar, br) {

	const dx = ax - bx, dy = ay - by, dz = az - bz;
	return Math.sqrt(dx*dx+dy*dy+dz*dz) - (ar+br);

};

xSSPS.prototype.updateRender = function(dt) {

	const hash3D = (x, y, z) => {
		const ix = (Math.floor(x/this.hSize)) + 1e5,
    	  	  iy = (Math.floor(y/this.hSize)) + 1e5,
    	  	  iz = (Math.floor(z/this.hSize)) + 1e5;
  		return ((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6;
	};

	const mHash3D = (x, y, z) => {
		const ix = (Math.floor(x/this.mhSize)) + 1e5,
    	  	  iy = (Math.floor(y/this.mhSize)) + 1e5,
    	  	  iz = (Math.floor(z/this.mhSize)) + 1e5;
  		return ((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6;
	};

	const hash = new Map();
	const mhash = new Map();

	const tkeys = [ 'lr', 'ldamp', 'gscale', 'opacity', 'incompress' ];

	// update state based on heat
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		const T = P.type;
		let found = false;
		for (let j=0; j<T.states.length; j++) {
			P.mesh[j].___wasVisible = !!P.mesh[j].___visible;
			P.mesh[j].___visible = false;
		}
		for (let j=0; j<T.states.length; j++) {
			const S = T.states[j];
			if (P.temp < S.ttemp1) {
				const mesh = P.mesh[j];
				for (let k=0; k<tkeys.length; k++) {
					const K = tkeys[k];
					P[K] = S[K];
				}
				mesh.material.opacity = P.opacity;
				mesh.___visible = true;
				break;
			}
			else if (P.temp < S.ttemp2 && j < (T.states.length-1)) {
				const t = Math.pow((P.temp - S.ttemp1) / (S.ttemp2 - S.ttemp1), 5.0);
				const S2 = T.states[j+1];
				const mesh = P.mesh[j];
				const mesh2 = P.mesh[j+1];
				for (let k=0; k<tkeys.length; k++) {
					const K = tkeys[k];
					P[K] = S[K] * (1-t) + S2[K] * t;
				}
				mesh.material.opacity = S.opacity * (1-t);
				mesh.___visible = true;
				if (t > 0.05) {
					mesh2.material.opacity = S2.opacity * t;
					mesh2.___visible = true;
				}
				break;
			}
		}
		for (let j=0; j<T.states.length; j++) {
			const mesh = P.mesh[j];
			if (mesh.___visible !== mesh.___wasVisible) {
				if (mesh.___visible) {
					this.scene.add(mesh);
				}
				else {
					this.scene.remove(mesh);
				}
			}
		}
		P.ldamp *= 0.25;
		window.AAAA = true;
	}

	// update position & angle, mark in hash
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		P.pos.x += P.vel.x * dt;
		P.pos.y += P.vel.y * dt;
		P.pos.z += P.vel.z * dt;
		P.ang.x += P.avel.x * dt;
		P.ang.y += P.avel.y * dt;
		P.ang.z += P.avel.z * dt;
		P.hkey = hash3D(P.pos.x, P.pos.y, P.pos.z);
		P.mhkey = mHash3D(P.pos.x, P.pos.y, P.pos.z);

		// insert into collision hash
		if (!hash.has(P.hkey)) {
			hash.set(P.hkey, [P]);
		}
		else {
			hash.get(P.hkey).push(P);
		}

		// insert into gravity hash
		if (!mhash.has(P.mhkey)) {
			mhash.set(P.mhkey, [P]);
		}
		else {
			mhash.get(P.mhkey).push(P);
		}

		// update THREE.js mesh
		for (let j=0; j<P.mesh.length; j++) {
			const mesh = P.mesh[j];
			if (mesh.___visible) {
				mesh.position.set(P.pos.x, P.pos.y, P.pos.z);
				if (j < 2) {
					mesh.rotation.set(P.ang.x, P.ang.y, P.ang.z, 'XYZ');
				}
				else {
					mesh.quaternion.copy(this.camera.quaternion);
				}
				const scale = P.lr;// * (P.gscale || 1);
				mesh.scale.set(scale, scale, scale);
				mesh.updateMatrix(true);
			}
		}
	}

	// compute gravity hash
	const gravity = [];
	for (let [hk, list] of mhash) {
		const gp = new THREE.Vector3(0, 0, 0);
		let tmass = 0.0;
		for (let i=0; i<list.length; i++) {
			const P = list[i];
			const mass = P.mass;
			gp.x += P.pos.x * mass;
			gp.y += P.pos.y * mass;
			gp.z += P.pos.z * mass;
			tmass += mass;
		}
		if (tmass > this.minMassGravity) {
			gp.x /= tmass;
			gp.y /= tmass;
			gp.z /= tmass;
			gravity.push({
				p: gp,
				mass: tmass
			});
		}
	}
	gravity.sort((a, b) => (b.mass - a.mass));
	gravity.length = Math.min(gravity.length, this.nGravity);

	let gtmass = gravity.length ? gravity[0].mass : 0;
	let gtp = gravity.length ? new THREE.Vector3(gravity[0].p.x, gravity[0].p.y, gravity[0].p.z) : new THREE.Vector3(0, 0, 0);
	for (let i=1; i<gravity.length; i++) {
		gtp.x = (gtp.x * gtmass + gravity[i].p.x * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtp.y = (gtp.y * gtmass + gravity[i].p.y * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtp.z = (gtp.z * gtmass + gravity[i].p.z * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtmass += gravity[i].mass;
	}

	if (!this.gtp) {
		this.gtp = gtp;
		this.gtmass = gtmass;
	}
	else {
		this.gtp.x += (gtp.x - this.gtp.x) * dt * 10;
		this.gtp.y += (gtp.y - this.gtp.y) * dt * 10;
		this.gtp.z += (gtp.z - this.gtp.z) * dt * 10;
		this.gtmass += (gtmass - this.gtmass) * dt * 6;
	}

	// update velocities with gravity
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		for (let j=0; j<gravity.length; j++) {
			const G = gravity[j];
			const dx = G.p.x - P.pos.x,
				  dy = G.p.y - P.pos.y,
				  dz = G.p.z - P.pos.z;
			const dlenSq = (dx*dx + dy*dy + dz*dz);
			if (dlenSq > 0.01) {
				const dlen1 = Math.sqrt(dlenSq);
				const dlen2 = G.mass / (Math.max(dlenSq, 20) * 0.1);
				P.vel.x += dlen2 * (dx / dlen1) * dt * this.gConst;
				P.vel.y += dlen2 * (dy / dlen1) * dt * this.gConst;
				P.vel.z += dlen2 * (dz / dlen1) * dt * this.gConst;
			}
		}
	}

	// handle collisions & sph
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		let ix, iy, iz;
		let dx, dy, dz, dlen;
		const tested = new Map();

		const fieldAdded = new Map();
		const fieldVel = new THREE.Vector3(0, 0, 0);
		let fieldWeight = 0;

		for (ix=-2; ix<=2; ix++) {
			for (iy=-2; iy<=2; iy++) {
				for (iz=-2; iz<=2; iz++) {
					const ihk = hash3D(
						P.pos.x + ix * this.hSize,
						P.pos.y + iy * this.hSize,
						P.pos.z + iz * this.hSize
					);
					if (hash.has(ihk)) {
						const list2 = hash.get(ihk);

						// calc field velocity
						for (let j=0; j<list2.length; j++) {
							const jP = list2[j];
							if (!fieldAdded.has(jP.id)) {
								fieldAdded.set(jP.id, true);
								const sdist = sphereDist(
									P.pos.x, P.pos.y, P.pos.z,
									jP.pos.x, jP.pos.y, jP.pos.z,
									this.hSize * 2, jP.lr
								);
								if (sdist < 0) {
									const tweight = jP.mass * Math.pow(1 - Math.max(-sdist / (this.hSize * 2), 1), 0.5);
									fieldWeight += tweight;
									fieldVel.x += jP.vel.x * tweight;
									fieldVel.y += jP.vel.y * tweight;
									fieldVel.z += jP.vel.z * tweight;
								}
							}
						}

						// hard collisions
						for (let j=0; j<list2.length; j++) {
							const jP = list2[j];
							if (jP.id > P.id && !tested.has(jP.id)) {
								tested.set(jP.id, true);
								const pR = P.lr * P.incompress * 0.75;
								const jpR = jP.lr * jP.incompress * 0.75;
								const sdist = sphereDist(
									P.pos.x, P.pos.y, P.pos.z,
									jP.pos.x, jP.pos.y, jP.pos.z,
									pR, jpR
								);
								if (sdist < 0) {
									dx = P.pos.x - jP.pos.x;
									dy = P.pos.y - jP.pos.y;
									dz = P.pos.z - jP.pos.z;
									dlen = Math.sqrt(dx*dx+dy*dy+dz*dz);
									if (dlen > 0.0001) {
										dx /= dlen;
										dy /= dlen;
										dz /= dlen;

										let x1;

										x1 = dx*P.vel.x + dy*P.vel.y + dz*P.vel.z;
										const vx1X = x1 * dx,
											  vx1Y = x1 * dy,
											  vx1Z = x1 * dz;
										const vy1X = P.vel.x - vx1X,
										      vy1Y = P.vel.y - vx1Y,
										      vy1Z = P.vel.z - vx1Z;
										const m1 = P.mass;

										x1 = -dx*P.vel.x + -dy*P.vel.y + -dz*P.vel.z;
										const vx2X = x1 * -dx,
											  vx2Y = x1 * -dy,
											  vx2Z = x1 * -dz;
										const vy2X = jP.vel.x - vx2X,
										      vy2Y = jP.vel.y - vx2Y,
										      vy2Z = jP.vel.z - vx2Z;
										const m2 = jP.mass;

										const TM = (m1+m2);
										const MA = (m1-m2)/TM,
											  MB = (m2-m1)/TM,
											  MC = (2*m2)/TM,
											  MD = (2*m1)/TM;

										const damp1 = Math.pow(P.ldamp, dt),
										      damp2 = Math.pow(jP.ldamp, dt);

										const ivx = (P.vel.x * m1 + jP.vel.x * m2) / TM,
											  ivy = (P.vel.y * m1 + jP.vel.y * m2) / TM,
											  ivz = (P.vel.z * m1 + jP.vel.z * m2) / TM;

										P.vel.x = (vx1X*MA + vx2X*MC + vy1X) * damp1 + ivx * (1-damp1);
										P.vel.y = (vx1Y*MA + vx2Y*MC + vy1Y) * damp1 + ivx * (1-damp1);
										P.vel.z = (vx1Z*MA + vx2Z*MC + vy1Z) * damp1 + ivx * (1-damp1);
										jP.vel.x = (vx1X*MD + vx2X*MB + vy2X) * damp2 + ivy * (1-damp2);
										jP.vel.y = (vx1Y*MD + vx2Y*MB + vy2Y) * damp2 + ivy * (1-damp2);
										jP.vel.z = (vx1Z*MD + vx2Z*MB + vy2Z) * damp2 + ivy * (1-damp2);
									}
								}
							}
						}
					}
				}
			}
		}

		// Apply field force
		if (fieldWeight > 0) {
			P.vel.x += ((P.vel.x * P.incompress + (fieldVel.x / fieldWeight) * (1 - P.incompress)) - P.vel.x) * dt * 8;
			P.vel.y += ((P.vel.y * P.incompress + (fieldVel.y / fieldWeight) * (1 - P.incompress)) - P.vel.y) * dt * 8;
			P.vel.z += ((P.vel.z * P.incompress + (fieldVel.z / fieldWeight) * (1 - P.incompress)) - P.vel.z) * dt * 8;
		}
	}

	return [{p: this.gtp, mass: this.gtmass}];

};

xSSPS.prototype.add = function(inArgs) {

	const args = inArgs || {};

	const rat = Math.PI * Math._random(), rap = Math.PI * 2.0 * Math._random();
	const rr = Math._random() * 100;
	const id = 1e9 + Math.floor(Math._random() * 1e9);
	const pos = args.pos ? args.pos : new THREE.Vector3(Math.cos(rat) * Math.sin(rap) * rr, Math.sin(rat) * Math.sin(rap) * rr, Math.cos(rap) * rr);

	let type = args.type || null;
	if (!type) {
		let tRandWeight = 0;
		for (let i=0; i<this.types.length; i++) {
			tRandWeight += this.types[i].randWeight;
		}
		const rand = Math._random() * tRandWeight;
		tRandWeight = 0;
		for (let i=0; i<this.types.length; i++) {
			if ((rand >= tRandWeight && rand < (tRandWeight + this.types[i].randWeight)) || (i===(this.types.length-1))) {
				type = this.types[i];
				break;
			}
			tRandWeight += this.types[i].randWeight;
		}
	}
	type = {...type, ...(type.__fn())};

	const temp = args.temp || 0;

	const ang = args.ang || new THREE.Vector3(Math._random()*50, Math._random()*50, Math._random()*50);

	const meshList = [];
	for (let i=0; i<type.states.length; i++) {
		const S = type.states[i];
		const mesh = new THREE.Mesh(S.geom, S.mat);
		mesh.position.set(pos.x, pos.y, pos.z);
		mesh.rotation.set(ang.x, ang.y, ang.z, 'XYZ');
		meshList.push(mesh);
	}

	const obj = {
		id: id,
		group: null,
		pos: pos,
		ang: ang,
		temp: 300,//temp,
		vel: args.vel || new THREE.Vector3(Math._random()*1-0.5, Math._random()*1-0.5, Math._random()*1-0.5),
		avel: args.avel || new THREE.Vector3(0, 0, 0),
		mass: type.mass,
		mesh: meshList,
		hr: 0.05,
		lr: 0.25,
		hdamp: 0.1,
		ldamp: 0.55,
		hkey: 0,
		type: type
	};

	this.list.push(obj);

};

xSSPS.prototype.seedTypes = function() {

	this.types = [];

	this.types.push({
		name: 'Hydrogen', id: 1,
		mass: 1,
		randWeight: 10,
		heatPerPressure: 1,
		heatDamp: 0.9,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999,
			lr: 0.25
		}, {
			ttemp1: 50,
			ttemp2: 60,
			lr: 0.5,
			ldamp: 0.25,
			clr: { r: 0, g: 0, b: 0.2 },
			eclr: { r: 0, g: 0, b: 0.05 },
			opacity: 0.5
		}, {
			ttemp1: 6900,
			ttemp2: 7000,
			lr: 0.4,
			ldamp: 0.5,
			clr: { r: 0.15, g: 0.15, b: 0.2 },
			opacity: 0.1
		}, {
			ttemp1: 7000,
			ttemp2: 1e10,
			lr: 3.5,
			ldamp: 0.85,
			clr: { r: 1.0, g: 0.5, b: 0.5 },
			opacity: 0.3
		}))
	});

	this.types.push({
		name: 'Iron', id: 2,
		mass: 15.0,
		randWeight: 2,
		heatPerPressure: 10,
		heatDamp: 0.1,
		__fn: () => (this.makeMaterial({
			ttemp1: 550,
			ttemp2: 600,
			ldamp: 0.05,
			lr: 0.6,
			ldamp: 0.85,
			clr: { r: 0.175, g: 0.175, b: 0.2 },
			eclr: { r: 0.175*.2, g: 0.175*.2, b: 0.2*.2 },
			opacity: 1.0
		}, {
			ttemp1: 600,
			ttemp2: 1e10,
			ldamp: 0.15,
			lr: 0.45,
			clr: { r: 0.85, g: 0.80, b: 0.0 },
			eclr: { r: 0.75, g: 0.70, b: 0.0 },
			opacity: 0.7
		}, {
			lr: 0.35,
			ttemp1: 1e10,
			ttemp2: 1e10 + 1e5
		}, {
			lr: 0.25,
			ttemp1: 1e11,
			ttemp2: 1e11 + 1e6
		}))
	});

	this.types.push({
		name: 'Water', id: 3,
		mass: 1,
		randWeight: 0.5,
		heatPerPressure: 1,
		heatDamp: 0.7,
		__fn: () => (this.makeMaterial({
			ttemp1: 300,
			ttemp2: 310,
			ldamp: 0.075,
			lr: 0.35,
			clr: { r: 0.3, g: 0.5, b: 0.8 },
			eclr: { r: 0.3/6, g: 0.5/6, b: 0.8/6 },
			opacity: 0.9
		}, {
			ttemp1: 400,
			ttemp2: 410,
			ldamp: 0.55,
			lr: 0.30,
			clr: { r: 0.3/6, g: 0.5/6, b: 0.3 },
			eclr: { r: 0.3/36, g: 0.5/36, b: 0.3/6 },
			opacity: 0.6
		}, {
			ttemp1: 1e10,
			ttemp2: 1e10 + 1e5,
			ldamp: 0.8,
			lr: 0.25,
			clr: { r: 0.8, g: 0.85, b: 1.0 },
			eclr: { r: 0.8/6, g: 0.85/6, b: 1.0/6 },
			opacity: 0.2
		}, {
			lr: 0.2,
			ttemp1: 1e11,
			ttemp2: 1e11 + 1e6
		}))
	});

	this.types.push({
		name: 'Nitrogen', id: 4,
		mass: 2,
		randWeight: 1,
		heatPerPressure: 1,
		heatDamp: 0.2,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999
		}, {
			ttemp1: 50,
			ttemp2: 60,
			lr: 0.5,
			ldamp: 0.25,
			clr: { r: 0, g: 0, b: 0.4 },
			eclr: { r: 0, g: 0, b: 0.1 },
			opacity: 0.75
		}, {
			ttemp1: 250,
			ttemp2: 1e10,
			lr: 0.75,
			ldamp: 0.5,
			clr: { r: 0.0, g: 0.0, b: 0.4 },
			opacity: 0.2
		}, {
			ttemp1: 1e10,
			ttemp2: 1e10 + 5,
			lr: 3.5,
			ldamp: 0.85,
			clr: { r: 1.0, g: 0.5, b: 1.0 },
			opacity: 0.3
		}))
	});

	this.types.push({
		name: 'Oxygen', id: 5,
		mass: 0.1,
		randWeight: 5,
		heatPerPressure: 5,
		heatDamp: 0.2,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999
		}, {
			ttemp1: 0,
			ttemp2: 10,
			lr: 0.35,
			ldamp: 0.25,
			clr: { r: 0.4, g: 0.4, b: 0.4 },
			eclr: { r: 0.4, g: 0.4, b: 0.4 },
			opacity: 0.75
		}, {
			ttemp1: 10,
			ttemp2: 5000,
			lr: 0.25,
			ldamp: 0.5,
			clr: { r: 0.7, g: 0.7, b: 0.7 },
			opacity: 0.2
		}, {
			ttemp1: 5000,
			ttemp2: 1e10 + 5,
			lr: 1.0,
			ldamp: 0.85,
			clr: { r: 1.0, g: 0.5, b: 1.0 },
			opacity: 0.3
		}))
	});

};

xSSPS.prototype.makeMaterial = function(iSolid, iLiquid, iGas, iPlasma) {

	const mkShaderGasMaterial = (r,g,b,a) => {
		const ret = new THREE.ShaderMaterial({
			uniforms: {
				color: { value: new THREE.Vector4(r, g, b, a) },
			},
			vertexShader: document.getElementById( 'vertexShader1' ).textContent,
			fragmentShader: document.getElementById( 'fragmentShader1' ).textContent,
			transparent: true
		});
		ret.depthTest = true;
		ret.depthWrite = false;
		ret.blending = THREE.AdditiveBlending;
		ret.needsUpdate = true;
		return ret;
	};

	const mkShaderSolidMaterial = (r,g,b,er,eg,eb,opacity) => {
		const ret = new THREE.MeshStandardMaterial({
			emissive: new THREE.Color(
				er >= 0 ? er : r*0.1,
				eg >= 0 ? eg : g*0.1,
				eb >= 0 ? eb : b*0.1
			),
			color: new THREE.Color(
				r, g, b
			),
			transparent: true
		});
		ret.opacity = opacity >= 0 ? opacity : 1.0;
		ret.depthTest = true;
		ret.depthWrite = false;
		ret.needsUpdate = true;
		return ret;
	};

	const solid = iSolid || {};
	const liquid = iLiquid || {};
	const gas = iGas || {};
	const plasma = iPlasma || {};

	return {
		states: [{
			ttemp1: solid.ttemp1 !== undefined ? solid.ttemp1 : 600,
			ttemp2: solid.ttemp2 !== undefined ? solid.ttemp2 : 650,
			lr: solid.lr || 0.5,
			ldamp: solid.ldamp || 0.05,
			mat: mkShaderSolidMaterial(
				solid.clr ? solid.clr.r : 0.2,
				solid.clr ? solid.clr.g : 0.2,
				solid.clr ? solid.clr.b : 0.2,
				solid.eclr ? solid.eclr.r : 0.02,
				solid.eclr ? solid.eclr.g : 0.02,
				solid.eclr ? solid.eclr.b : 0.02,
				solid.opacity >= 0 ? solid.opacity : 0.9
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 20/60,
			opacity: solid.opacity >= 0 ? solid.opacity : 0.9,
			geom: this.pgeom2,
			gscale: 1.0
		},
		{
			ttemp1: liquid.ttemp1 !== undefined ? liquid.ttemp1 : 1200,
			ttemp2: liquid.ttemp2 !== undefined ? liquid.ttemp2 : 1300,
			lr: liquid.lr || 0.65,
			ldamp: liquid.ldamp || 0.35,
			mat: mkShaderSolidMaterial(
				liquid.clr ? liquid.clr.r : 0.5,
				liquid.clr ? liquid.clr.g : 0.2,
				liquid.clr ? liquid.clr.b : 0.0,
				liquid.eclr ? liquid.eclr.r : 0.02,
				liquid.eclr ? liquid.eclr.g : 0.02,
				liquid.eclr ? liquid.eclr.b : 0.02,
				liquid.opacity >= 0 ? liquid.opacity : 0.4
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 5/60,
			opacity: liquid.opacity >= 0 ? liquid.opacity : 0.4,
			geom: this.pgeom2,
			gscale: 1.2
		},
		{
			ttemp1: gas.ttemp1 !== undefined ? gas.ttemp1 : 3300,
			ttemp2: gas.ttemp2 !== undefined ? gas.ttemp2 : 3350,
			lr: gas.lr || 1.5,
			ldamp: gas.ldamp || 0.75,
			mat: mkShaderGasMaterial(
				gas.clr ? gas.clr.r : 1.0,
				gas.clr ? gas.clr.g : 0.5,
				gas.clr ? gas.clr.b : 0.05,
				gas.opacity >= 0 ? gas.opacity : 0.3
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 1/60,
			opacity: gas.opacity >= 0 ? gas.opacity : 0.3,
			geom: this.pgeom,
			gscale: 1.3
		},
		{
			ttemp1: plasma.ttemp1 !== undefined ? plasma.ttemp1 : 8300,
			ttemp2: plasma.ttemp2 !== undefined ? plasma.ttemp2 : 8800,
			lr: plasma.lr || 4.5,
			ldamp: plasma.ldamp || 0.85,
			mat: mkShaderGasMaterial(
				plasma.clr ? plasma.clr.r : 1.0,
				plasma.clr ? plasma.clr.g : 1.0,
				plasma.clr ? plasma.clr.b : 1.0,
				plasma.opacity >= 0 ? plasma.opacity : 0.3,
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 30/60,
			opacity: plasma.opacity >= 0 ? plasma.opacity : 0.3,
			geom: this.pgeom,
			gscale: 1.4
		}]
	};
};