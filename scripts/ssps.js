window.xSSPS = function(PCOUNT, RSIZE) {

    this.RSIZE = RSIZE || 512;
    this.PCOUNT = PCOUNT || 512;

    this.hashLength = 128;
    this.hashWSize = 3.0;
    this.hashBinLen = 128;
    this.rayCount = 128;
    this.rayDist = 20;
    this.rfScale = 2.0;

    this.restDensity = 0.1;
    this.fieldLen = 4;
    this.gConst = 0.02;
    this.bound = 15;

    this.shootIndex = 0;
    this.lkeys = {};

    this.canvas = document.createElement('canvas');
    this.ctx3d = this.canvas.getContext('webgl');
    this.gpu = new GPU({
        canvas: this.canvas,
        webGl: this.ctx3d,
        mode: 'webgl'
    });

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
        
        up = normalize(cross(cross(cd, up), cd));

        vec3 vup = normalize(cross(cd, up));
        vec3 vleft = up;

        vec2 X = (uv - vec2(0.5, 0.5)) * s0;

        return c0 + vup * X.y + vleft * X.x;

    }`);

    this.gpu.addNativeFunction('getRayDir', `vec3 getRayDir(vec3 R0,  vec2 uv,  vec3 c0, vec3 cd, vec3 up, float s0, float s1, float slen) {
        
        up = normalize(cross(cross(cd, up), cd));

        vec3 vup = normalize(cross(cd, up));
        vec3 vleft = up;

        vec2 X = (uv - vec2(0.5, 0.5)) * s1;

        vec3 far = vup * X.y + vleft * X.x;

        return normalize((normalize(cd) * slen + far + c0) - R0);

    }`);

    this.gpu.addNativeFunction('reflectWrap', `vec3 reflectWrap(vec3 I, vec3 N) {

        return normalize(reflect(normalize(I), normalize(N)));
    
    }`);

    this.gpu.addNativeFunction('refractWrap', `vec3 refractWrap(vec3 I, vec3 N, float eta) {

        return normalize(refract(normalize(I), normalize(N), eta));
    
    }`);

    this.gpu.addNativeFunction('distanceWrap', `float distanceWrap(vec3 A, vec3 B) {

        return length(A - B);
    
    }`);

    /*
        UPDATE SUPPORT FUNCTIONS (GLSL)
     */

    this.gpu.addNativeFunction('compGravity', `vec3 compGravity(vec3 M, vec3 J, float jMass, float f) {

        vec3 delta = J - M;
        float dlen = length(delta);
                 
        if (dlen > 0.05) {
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

    this.velocityKernel = this.gpu.createKernel(function(positions, velocities, attrs, pullMass, playerPos, oDensity, oVisc, oMass, oIncomp, dt) {

        // Unpack particle
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
        var mdensity = oDensity;
        var mmass = oMass;
        var incomp = oIncomp;
        var viscdt = oVisc;
        var ddensity = 0.0,
            nddensity = 0.0;

        // Compute pressure on this particle
        for (var i=0; i<this.constants.PCOUNT; i++) {
            var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
            var density = oDensity;
            var dret = [0, 0]; dret = partDensity(mpos, opos, this.constants.fieldLen);
            ddensity += dret[0] * density;
            nddensity += dret[1] * density;
        }

        // Interact with player by adding pressure
        var opos = [playerPos[0], playerPos[1], playerPos[2]];
        var fl2 = this.constants.fieldLen;
        var pf = [0, 0];
        if (pullMass > 0.0) {
            pf = partDensity(mpos, opos, fl2);
        }
        ddensity += pf[0] * (pullMass / 10.0 * this.constants.restDensity);
        nddensity += pf[1] * (pullMass / 10.0 * this.constants.restDensity);

        var spressure = (ddensity - this.constants.restDensity) * incomp;
        var snpressure = nddensity * incomp;

        // Compute force from pressure on this particle
        var ret = [mvx, mvy, mvz];

        for (var i=0; i<this.constants.PCOUNT; i++) {
            if (Math.abs(i - me) > 0.01) {
                var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
                var ovel = [velocities[i*3], velocities[i*3+1], velocities[i*3+2]];
                var mass = oMass;
                var jmf = (2. * mass) / (mass + mmass);
                var dret = [0., 0., 0.]; dret = pressForce(mpos, opos, mvel, ovel, this.constants.fieldLen, dt, spressure, snpressure, viscdt);
                ret[0] += dret[0] * jmf; ret[1] += dret[1] * jmf; ret[2] += dret[2] * jmf;
            }
        }

        // Player pressure
        if (pullMass > 0.0) {
            var opos = [playerPos[0], playerPos[1], playerPos[2]];
            var ovel = [0, 0, 0];
            var mass = 1.0;
            var jmf = (2. * mass) / (mass + mmass);
            var dret = [0., 0., 0.]; dret = pressForce(mpos, opos, mvel, ovel, this.constants.fieldLen, dt, spressure, snpressure, viscdt);
            ret[0] += dret[0] * jmf; ret[1] += dret[1] * jmf; ret[2] += dret[2] * jmf;          
        }

        // Compute gravitational force on this particle
        var gf = this.constants.gConst * dt * 0.1;
        for (var i=0; i<this.constants.PCOUNT; i++) {
            var opos = [positions[i*3], positions[i*3+1], positions[i*3+2]];
            var mass = attrs[i*6+2];
            var dret = [0., 0., 0.]; dret = compGravity(mpos, opos, mass, gf);
            ret[0] += dret[0]; ret[1] += dret[1]; ret[2] += dret[2];
        }

        // Enforce scene boundaries
        if (ret[0] < 0 && mpos[0] < -this.constants.bound) {
            ret[0] = -ret[0];
        }
        else if (ret[0] > 0 && mpos[0] > this.constants.bound) {
            ret[0] = -ret[0];
        }

        if (ret[1] < 0 && mpos[1] < -this.constants.bound) {
            ret[1] = -ret[1];
        }
        else if (ret[1] > 0 && mpos[1] > this.constants.bound) {
            ret[1] = -ret[1];
        }

        if (ret[2] < 0 && mpos[2] < -this.constants.bound) {
            ret[2] = -ret[2];
        }
        else if (ret[2] > 0 && mpos[2] > this.constants.bound) {
            ret[2] = -ret[2];
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
            bound: this.bound,
            fieldLen: this.fieldLen,
            gConst: this.gConst,
            restDensity: ((4 / 3) * Math.PI * Math.pow(this.fieldLen, 3)) * this.restDensity
        },
        loopMaxIterations: this.PCOUNT,
        output: [ this.PCOUNT*3 ],
        outputToTexture: true,
        paramTypes: {
            particles: 'NumberTexture',
            velocities: 'NumberTexture',
            attrs: 'NumberTexture',
            playerPos: 'Array(3)',
            pullMass: 'Number',
            dt: 'Number'
        }
    });

    /*
        SET POSITION/VELOCITES
     */

    const mkSetter = (pitch) => {
        return this.gpu.createKernel(function(list, index, setv) {
            if (Math.abs(index - Math.floor(this.thread.x / this.constants.pitch)) < 0.01) {
                return setv[this.thread.x % this.constants.pitch];
            }
            else {
                return list[this.thread.x];
            }
        }, {
            output: [ this.PCOUNT*pitch ],
            outputToTexture: true,
            constants: { PCOUNT: this.PCOUNT, pitch },
            loopMaxIterations: this.PCOUNT,
            paramTypes: {
                list: 'NumberTexture',
                index: 'Number',
                setv: 'Array(3)'
            }
        });
    };
    
    this.setPosKernel = mkSetter(3);
    this.setVelKernel = mkSetter(3);
    this.setAttrKernel = mkSetter(6);

    /*
        POSITION UPDATE KERNEL
     */

    this.positionKernel = this.gpu.createKernel(function(positions, velocities, dt) {
        return positions[this.thread.x] + velocities[this.thread.x] * dt;
    }, {
        output: [ this.PCOUNT*3 ],
        outputToTexture: true,
        loopMaxIterations: this.PCOUNT,
        paramTypes: {
            particles: 'NumberTexture',
            velocities: 'NumberTexture',
            dt: 'Number'
        }
    });

    this.camRotKernel = this.gpu.createKernel(function(rotLR, rotUD, camCenter, camDir, camUp, camNearWidth, camFarWidth, camDist) {
        var uv = [0.5+rotLR, 0.5+rotUD];

        var CC = [camCenter[0], camCenter[1], camCenter[2]],
            CD = [camDir[0], camDir[1], camDir[2]],
            CU = [camUp[0], camUp[1], camUp[2]];

        var ray0 = [0,0,0]; ray0 = getRay0(uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);
        var rayDir = [0,0,0]; rayDir = getRayDir(ray0, uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);

        if (this.thread.x < 0.01) {
            return rayDir[0];
        }
        else if (this.thread.x < 1.01) {
            return rayDir[1];
        }
        else if (this.thread.x < 2.01) {
            return rayDir[2];
        }
    }, {
        output: [ 3 ],
        paramTypes: {
            rotLR: 'Number',
            rotUD: 'Number',
            camCenter: 'Array(3)',
            camDir: 'Array(3)',
            camUp: 'Array(3)',
            camNearWidth: 'Number',
            camFarWidth: 'Number',
            camDist: 'Number'
        }
    });

    /*
        RENDER KERNEL
     */
    
    this.renderMode = 0;
    
    this.renderKernel = [
        this.gpu.createKernel(function(hash, camCenter, camDir, camUp, camNearWidth, camFarWidth, camDist) {

            var uv = [this.thread.x / (this.constants.RSIZE-1), this.thread.y / (this.constants.RSIZE-1)];

            var CC = [camCenter[0], camCenter[1], camCenter[2]],
                CD = [camDir[0], camDir[1], camDir[2]],
                CU = [camUp[0], camUp[1], camUp[2]];

            var ray0 = [0,0,0]; ray0 = getRay0(uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);
            var rayDir = [0,0,0]; rayDir = getRayDir(ray0, uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);

            var oray0 = [ray0[0], ray0[1], ray0[2]];

            this.color(0., 0., 0., 1.);

            var outClr = [
                0., 0., 0.
            ];

            var int = 0., light = 0.;

            var stepR = this.constants.rayDist / this.constants.rayCount;
            var max2 = Math.floor(this.constants.bound * 2.);
            var found = -1.0;
            var norm = [0., 0., 0.];
            ray0[0] += rayDir[0] * stepR * 10.;
            ray0[1] += rayDir[1] * stepR * 10.;
            ray0[2] += rayDir[2] * stepR * 10.;
            for (var i=0; i<this.constants.rayCount; i++) {
                if (found < 0.) {
                    ray0[0] += rayDir[0] * stepR;
                    ray0[1] += rayDir[1] * stepR;
                    ray0[2] += rayDir[2] * stepR;

                    var ix = (Math.floor(ray0[0]/this.constants.hashWSize)) + max2,
                        iy = (Math.floor(ray0[1]/this.constants.hashWSize)) + max2, 
                        iz = (Math.floor(ray0[2]/this.constants.hashWSize)) + max2; 
                    var k = (Math.round(ix * 137.) + Math.round(iy * 197.) + Math.round(iz * 167.)) % this.constants.hashLength;

                    var fq = 0.0;

                    norm[0] = 0.;
                    norm[1] = 0.;
                    norm[2] = 0.;

                    var off = Math.round(k * this.constants.hashBinLen * 4.);
                    var endj = -1.0;
                    for (var j=0; j<this.constants.hashBinLen; j++) {
                        if (endj < 0.) {
                            if (hash[off+3.] < 0.5) {
                                endj = 1.0;
                            }
                            else {
                                var dx = hash[off], dy = hash[off+1.], dz = hash[off+2.];
                                dx -= ray0[0]; dy -= ray0[1]; dz -= ray0[2];
                                var r = dx*dx + dy*dy + dz*dz;
                                if (r > 0.0) {
                                    dx /= r;
                                    dy /= r;
                                    dz /= r;

                                    r = Math.sqrt(r);
                                    r /= this.constants.rfScale;
                                    if (r <= 0.707) {
                                        var q = (r*r*r*r - r*r + 0.25);
                                        norm[0] += dx*q;
                                        norm[1] += dy*q;
                                        norm[2] += dz*q;
                                        fq += q;
                                    }
                                }
                                off += 4.;
                            }
                        }
                    }

                    if (fq > 0.175) {
                        norm[0] /= fq;
                        norm[1] /= fq;
                        norm[2] /= fq;
                        var ref = [0, 0, 0]; ref = refractWrap(norm, oray0, 1./1.5);
                        var dot = ref[0] * rayDir[0] + ref[1] * rayDir[1] + ref[2] * rayDir[2];
                        int = (1. - Math.pow(i / this.constants.rayCount, 3.0)) * 0.7;
                        light = Math.min(1., Math.max(dot*dot*dot*dot, 0.)) * int;
                        found = 1.0;
                    }
                }
            }

            outClr[2] = int + light;
            outClr[1] = light;
            outClr[0] = light;

            this.color(Math.min(outClr[0], 1.), Math.min(outClr[1], 1.), Math.min(outClr[2], 1.), 1.);

        }, {
            graphical: true,
            constants: {
                hashLength: this.hashLength,
                hashWSize: this.hashWSize,
                hashBinLen: this.hashBinLen,
                rayCount: this.rayCount,
                rayDist: this.rayDist,
                rfScale: this.rfScale,
                bound: this.bound,
                RSIZE: this.RSIZE,
                PCOUNT: this.PCOUNT,
            },
            loopMaxIterations: this.rayCount * this.hashBinLen,
            output: [this.RSIZE, this.RSIZE],
            paramTypes: {
                hash: 'Array',
                attrs: 'NumberTexture',
                camCenter: 'Array(3)',
                camDir: 'Array(3)',
                camUp: 'Array(3)',
                camNearWidth: 'Number',
                camFarWidth: 'Number',
                camDist: 'Number'
            }
        }),
        this.gpu.createKernel(function(hash, camCenter, camDir, camUp, camNearWidth, camFarWidth, camDist) {

            var uv = [this.thread.x / (this.constants.RSIZE-1), this.thread.y / (this.constants.RSIZE-1)];

            var CC = [camCenter[0], camCenter[1], camCenter[2]],
                CD = [camDir[0], camDir[1], camDir[2]],
                CU = [camUp[0], camUp[1], camUp[2]];

            var ray0 = [0,0,0]; ray0 = getRay0(uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);
            var rayDir = [0,0,0]; rayDir = getRayDir(ray0, uv, CC, CD, CU, camNearWidth, camFarWidth, camDist);

            var oray0 = [ray0[0], ray0[1], ray0[2]];

            this.color(0., 0., 0., 1.);

            var outClr = [
                0., 0., 0.
            ];

            var int = 0., light = 0.;

            var stepR = this.constants.rayDist / this.constants.rayCount;
            var max2 = Math.floor(this.constants.bound * 2.);
            var found = -1.0;
            var norm = [0., 0., 0.];
            ray0[0] += rayDir[0] * stepR * 10.;
            ray0[1] += rayDir[1] * stepR * 10.;
            ray0[2] += rayDir[2] * stepR * 10.;
            for (var i=0; i<this.constants.rayCount; i++) {
                if (found < 0.) {
                    ray0[0] += rayDir[0] * stepR;
                    ray0[1] += rayDir[1] * stepR;
                    ray0[2] += rayDir[2] * stepR;

                    var ix = (Math.floor(ray0[0]/this.constants.hashWSize)) + max2,
                        iy = (Math.floor(ray0[1]/this.constants.hashWSize)) + max2, 
                        iz = (Math.floor(ray0[2]/this.constants.hashWSize)) + max2; 
                    var k = (Math.round(ix * 137.) + Math.round(iy * 197.) + Math.round(iz * 167.)) % this.constants.hashLength;

                    var fq = 0.0;

                    norm[0] = 0.;
                    norm[1] = 0.;
                    norm[2] = 0.;

                    var off = Math.round(k * this.constants.hashBinLen * 4.);
                    var endj = -1.0;
                    for (var j=0; j<this.constants.hashBinLen; j++) {
                        if (endj < 0.) {
                            if (hash[off+3.] < 0.5) {
                                endj = 1.0;
                            }
                            else {
                                var dx = hash[off], dy = hash[off+1.], dz = hash[off+2.];
                                dx -= ray0[0]; dy -= ray0[1]; dz -= ray0[2];
                                var r = dx*dx + dy*dy + dz*dz;
                                if (r > 0.0) {
                                    dx /= r;
                                    dy /= r;
                                    dz /= r;

                                    r = Math.sqrt(r);
                                    r /= this.constants.rfScale;
                                    if (r <= 0.707) {
                                        var q = (r*r*r*r - r*r + 0.25);
                                        norm[0] += dx*q;
                                        norm[1] += dy*q;
                                        norm[2] += dz*q;
                                        fq += q;
                                    }
                                }
                                off += 4.;
                            }
                        }
                    }

                    if (fq > 0.175) {
                        norm[0] /= fq;
                        norm[1] /= fq;
                        norm[2] /= fq;
                        var ref = [0, 0, 0]; ref = refractWrap(norm, oray0, 1./1.5);
                        var dot = ref[0] * rayDir[0] + ref[1] * rayDir[1] + ref[2] * rayDir[2];
                        int = (1. - Math.pow(i / this.constants.rayCount, 3.0)) * 0.7;
                        light = Math.min(1., Math.max(dot*dot*dot*dot, 0.)) * int;
                        found = 1.0;
                    }
                }
            }

            outClr[2] = int + light;
            outClr[1] = light;
            outClr[0] = light;

            if (found > 0.) {
                found = -1.;
                var nrayDir = [-norm[0], -norm[1], -norm[2]];
                var nstepR = this.constants.rayDistRef / this.constants.rayCountRef;
                for (var i=0; i<this.constants.rayCountRef; i++) {
                    if (found < 0.) {
                        ray0[0] += nrayDir[0] * nstepR;
                        ray0[1] += nrayDir[1] * nstepR;
                        ray0[2] += nrayDir[2] * nstepR;

                        var ix = (Math.floor(ray0[0]/this.constants.hashWSize)) + max2,
                            iy = (Math.floor(ray0[1]/this.constants.hashWSize)) + max2, 
                            iz = (Math.floor(ray0[2]/this.constants.hashWSize)) + max2; 
                        var k = (Math.round(ix * 137.) + Math.round(iy * 197.) + Math.round(iz * 167.)) % this.constants.hashLength;

                        var fq = 0.0;

                        norm[0] = 0.;
                        norm[1] = 0.;
                        norm[2] = 0.;

                        var off = Math.round(k * this.constants.hashBinLen * 4.);
                        var endj = -1.0;
                        for (var j=0; j<this.constants.hashBinLen; j++) {
                            if (endj < 0.) {
                                if (hash[off+3.] < 0.5) {
                                    endj = 1.0;
                                }
                                else {
                                    var dx = hash[off], dy = hash[off+1.], dz = hash[off+2.];
                                    dx -= ray0[0]; dy -= ray0[1]; dz -= ray0[2];
                                    var r = dx*dx + dy*dy + dz*dz;
                                    if (r > 0.0) {
                                        dx /= r;
                                        dy /= r;
                                        dz /= r;

                                        r = Math.sqrt(r);
                                        r /= this.constants.rfScale;
                                        if (r <= 0.707) {
                                            var q = (r*r*r*r - r*r + 0.25);
                                            norm[0] += dx*q;
                                            norm[1] += dy*q;
                                            norm[2] += dz*q;
                                            fq += q;
                                        }
                                    }
                                    off += 4.;
                                }
                            }
                        }

                        if (fq > 0.175) {
                            norm[0] /= fq;
                            norm[1] /= fq;
                            norm[2] /= fq;
                            var ref = [0, 0, 0]; ref = refractWrap(norm, oray0, 1./1.5);
                            var dot = ref[0] * rayDir[0] + ref[1] * rayDir[1] + ref[2] * rayDir[2];
                            int = (1. - Math.pow(i / this.constants.rayCountRef, 3.0));
                            light = Math.min(1., Math.max(dot*dot*dot*dot, 0.)) * int;
                            found = 1.0;
                        }
                    }
                }

                if (found > 0.) {
                    outClr[2] += (int + light) * 0.1;
                    outClr[1] += (light) * 0.1;
                    outClr[0] += (light) * 0.1;
                }
            }

            this.color(Math.min(outClr[0], 1.), Math.min(outClr[1], 1.), Math.min(outClr[2], 1.), 1.);

        }, {
            graphical: true,
            constants: {
                hashLength: this.hashLength,
                hashWSize: this.hashWSize,
                hashBinLen: this.hashBinLen,
                rayCount: this.rayCount,
                rayDist: this.rayDist,
                rayCountRef: Math.floor(this.rayCount / 4),
                rayDistRef: Math.floor(this.rayDist / 4),
                rfScale: this.rfScale,
                bound: this.bound,
                RSIZE: this.RSIZE,
                PCOUNT: this.PCOUNT,
            },
            loopMaxIterations: this.rayCount * this.hashBinLen,
            output: [this.RSIZE, this.RSIZE],
            paramTypes: {
                hash: 'Array',
                attrs: 'NumberTexture',
                camCenter: 'Array(3)',
                camDir: 'Array(3)',
                camUp: 'Array(3)',
                camNearWidth: 'Number',
                camFarWidth: 'Number',
                camDist: 'Number'
            }
        })
    ];

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
            output: [size]
        });
        return kern(arr);
    };

    const seedPositions = this.gpu.createKernel(function(){
        var t = random(Math.floor(this.thread.x / 3), this.constants.seed + 0.5);
        if (Math.floor(this.thread.x / 3) > (this.constants.PCOUNT-33)) {
            t + 1.0;
        }
        var r = t * this.constants.maxr;
        var a = 3.141592 * 2.0 * random(Math.floor(this.thread.x / 3), this.constants.seed + 1.5);
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
            maxr: 25,
            seed: Math.random() * 1e6,
            PCOUNT: this.PCOUNT
        },
        outputToTexture: true,
        output: [this.PCOUNT * 3]
    });

    const seedVelocities = this.gpu.createKernel(function(){
        var t = random(this.thread.x, this.constants.seed + 0.5);
        return (t - 0.5) * this.constants.iv;
    }, {
        constants: {
            iv: 1.5,
            seed: Math.random() * 1e6,
            PCOUNT: this.PCOUNT
        },
        outputToTexture: true,
        output: [this.PCOUNT * 3]
    });

    this.sDensity = [0.125, 0.25, 0.5, 1.0, 1.5, 2.0]; this.sDensityI = 5;
    this.sMass = [0.05, 0.1, 0.2, 0.4, 0.8]; this.sMassI = 2;
    this.sIncomp = [1.0, 0.8, 0.6, 0.4, 0.2]; this.sIncompI = 3;
    this.sVisc = [0.05, 0.1, 0.25, 0.35, 0.5, 0.75, 0.95]; this.sViscI = 0;

    const seedAttrs = this.gpu.createKernel(function(){
        var comp = this.thread.x % 6;
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
        else if (comp < 4.01) {
            return 0.35; // visc
        }
        else {
            // type
            if (Math.floor(this.thread.x / 6) > (this.constants.PCOUNT-33)) {
                return 1.0;
            }
            else {
                return 0.0;
            }
        }
    }, {
        constants: {
            iv: 50,
            seed: Math.random() * 1e6,
            PCOUNT: this.PCOUNT
        },
        outputToTexture: true,
        output: [this.PCOUNT * 6]
    });

    /*
        COPPIER
     */
    
    const makeCoppier = (isize) => {
        return this.gpu.createKernel(function(parray) {
            return parray[this.thread.x];
        }, {
            outputToTexture: true,
            output: [isize * this.PCOUNT]
        })
    };

    /*
     * Seed Particles
     */ 
    this.reset = () => {
        this.cam = {
            p: {x: 0, y: 0, z: 17.5},
            dir: {x: 0, y: 0, z: -1},
            up: {x: 0, y: 1, z: 0},
            nearW: 0.0001,
            farW: 1000,
            farDist: 1000
        };

        this.move = {
            toLR: 0,
            toUD: 0,
            toFR: 0,
            tLR: 0,
            tUD: 0,
            tFR: 0,
            toPull: 0,
            tPull: 0,
            toPullR: 0,
            tPullR: 0
        };

        this.data = {
            pos: seedPositions(),
            vel: seedVelocities(),
            attr: seedAttrs()
        };
    };

    this.reset();

    /*
        Make coppiers
     */
    this.posCoppier = makeCoppier(3);
    this.velCoppier = makeCoppier(3);
    this.attrCoppier = makeCoppier(6);

    this.readArray = this.gpu.createKernel(function(arr){
        return arr[this.thread.x];
    }, {
        outputToTexture: false,
        output: [this.PCOUNT * 3]
    });

    this.hash = [];
    for (let i=0; i<this.hashLength; i++) {
        for (let j=0; j<this.hashBinLen; j++) {
            this.hash.push(0); this.hash.push(0); this.hash.push(0); this.hash.push(-1);
        }
    }

};

xSSPS.prototype.hashFn = function(x, y, z) {
    const max2 = Math.floor(this.bound * 2);
    const ix = (Math.floor(x/this.hashWSize)) + max2,
          iy = (Math.floor(y/this.hashWSize)) + max2, 
          iz = (Math.floor(z/this.hashWSize)) + max2; 
    return ((ix * 137) + (iy * 197) + (iz * 167)) % this.hashLength;
};

xSSPS.prototype.updateRender = function(keys, dt) {

    this.handleInput(keys, dt);

    this.sDensityI = this.sDensityI % this.sDensity.length;
    this.sViscI = this.sViscI % this.sVisc.length;
    this.sMassI = this.sMassI % this.sMass.length;
    this.sIncompI = this.sIncompI % this.sIncomp.length;

    const SUBSTEPS = 2;

    for (let i=0; i<SUBSTEPS; i++) {
        // Update velocities via gravity & pressure
        this.data.vel = this.velocityKernel(
            this.data.pos,
            this.velCoppier(this.data.vel),
            this.data.attr,
            this.move.tPull,
            [this.cam.p.x + this.cam.dir.x * this.move.tPullR, this.cam.p.y + this.cam.dir.y * this.move.tPullR, this.cam.p.z + this.cam.dir.z * this.move.tPullR],
            this.sDensity[this.sDensityI],
            this.sVisc[this.sViscI],
            this.sMass[this.sMassI],
            this.sIncomp[this.sIncompI],
            dt / SUBSTEPS
        );

        // Update positions
        this.data.pos = this.positionKernel(
            this.posCoppier(this.data.pos),
            this.data.vel,
            dt / SUBSTEPS
        );
    }

    // Update hash for raytracing
    // 

    // Clear hash
    for (let i=0; i<this.hashLength; i++) {
        const o1 = i * this.hashBinLen * 4;
        for (let j=0; j<this.hashBinLen; j++) {
            const o2 = o1 + j * 4;
            this.hash[o2 + 3] = 0.;
        }
    }

    // Insert points into hash
    const positions = this.readArray(this.data.pos);
    const hmap = {};
    for (let i=0; i<this.PCOUNT; i++) {
        const off = i * 3;
        const x = positions[off], y = positions[off+1], z = positions[off+2];
        const r = this.rfScale * 0.707;
        const ix1 = Math.floor((x - r) / this.hashWSize) - 1,
              iy1 = Math.floor((y - r) / this.hashWSize) - 1,
              iz1 = Math.floor((z - r) / this.hashWSize) - 1,
              ix2 = Math.round((x + r) / this.hashWSize) + 1,
              iy2 = Math.round((y + r) / this.hashWSize) + 1,
              iz2 = Math.round((z + r) / this.hashWSize) + 1;
        for (let ix=ix1; ix<=ix2; ix ++) {
            for (let iy=iy1; iy<=iy2; iy ++) {
                for (let iz=iz1; iz<=iz2; iz ++) {
                    const xc = (ix + 0.5) * this.hashWSize,
                          yc = (iy + 0.5) * this.hashWSize,
                          zc = (iz + 0.5) * this.hashWSize;
                    const dx = xc - x, dy = yc - y, dz = zc - z;
                    const dist = Math.sqrt(dx*dx+dy*dy+dz*dz);
                    if (dist <= (r+this.hashWSize)) {
                        const hkey = this.hashFn(xc, yc, zc);
                        (hmap[hkey] = hmap[hkey] || []).push({x, y, z, type: i < ((this.PCOUNT-1)-32) ? 1. : 2., dist})
                    }
                }
            }
        }
    }
    const sortFn = (a, b) => (a.dist - b.dist);
    for (var hkeyStr in hmap) {
        const hkey = parseInt(hkeyStr);
        const list = hmap[hkeyStr];
        list.sort(sortFn);
        for (let i=0; i<list.length && i<this.hashBinLen; i++) {
            const P = list[i];
            const off = hkey * this.hashBinLen * 4 + i * 4;
            this.hash[off+0] = P.x;
            this.hash[off+1] = P.y;
            this.hash[off+2] = P.z;
            this.hash[off+3] = P.type;
        }
    }

    // Render
    this.renderMode = this.renderMode % this.renderKernel.length;
    this.renderKernel[this.renderMode](
        this.hash,
        [this.cam.p.x, this.cam.p.y, this.cam.p.z],
        [this.cam.dir.x, this.cam.dir.y, this.cam.dir.z],
        [this.cam.up.x, this.cam.up.y, this.cam.up.z],
        this.cam.nearW, this.cam.farW, this.cam.farDist
    );

};

xSSPS.prototype.handleInput = function(keys, dt) {

    this.move.toLR = this.move.toFR = this.move.toUD = 0;

    if (keys[37]) {
        this.move.toLR -= 1;
    }
    if (keys[39]) {
        this.move.toLR += 1;
    }
    if (keys[38]) {
        this.move.toUD -= 1;
    }
    if (keys[40]) {
        this.move.toUD += 1;
    }
    if (keys[83]) {
        this.move.toFR -= 1;
    }
    if (keys[87]) {
        this.move.toFR += 1;
    }

    if (this.lkeys[82] && !keys[82]) {
        this.renderMode += 1;
    }
    if (this.lkeys[27] && !keys[27]) {
        this.reset();
    }
    if (this.lkeys[49] && !keys[49]) {
        this.sDensityI += 1;
    }
    if (this.lkeys[50] && !keys[50]) {
        this.sViscI += 1;
    }
    if (this.lkeys[51] && !keys[51]) {
        this.sMassI += 1;
    }
    if (this.lkeys[52] && !keys[52]) {
        this.sIncompI += 1;
    }

    if (this.lkeys[32] && !keys[32]) {
        this.data.pos = this.setPosKernel(
            this.posCoppier(this.data.pos),
            (this.PCOUNT-1) - this.shootIndex,
            [
                this.cam.p.x + this.cam.dir.x * 0.5,
                this.cam.p.y + this.cam.dir.y * 0.5,
                this.cam.p.z + this.cam.dir.z * 0.5
            ]
        );
        this.data.vel = this.setVelKernel(
            this.velCoppier(this.data.vel),
            (this.PCOUNT-1) - this.shootIndex,
            [
                this.cam.dir.x * 15,
                this.cam.dir.y * 15,
                this.cam.dir.z * 15
            ]
        );
        this.shootIndex = (this.shootIndex + 1) % 32;
    }
    if (keys[69]) {
        this.move.toPull = 10;
        this.move.toPullR = 5;
    }
    else {
        this.move.toPull = -1;
        this.move.toPullR = 0.25;
    }

    this.lkeys = {...keys};

    this.move.tLR += (this.move.toLR - this.move.tLR) * dt * 1.5;
    this.move.tUD += (this.move.toUD - this.move.tUD) * dt * 1.5;
    this.move.tFR += (this.move.toFR - this.move.tFR) * dt * 1.5;
    this.move.tPull += (this.move.toPull - this.move.tPull) * dt * 1.5;
    this.move.tPullR += (this.move.toPullR - this.move.tPullR) * dt * 1.5;

    const utmp = squat.normv(squat.crossv(squat.crossv([this.cam.dir.x, this.cam.dir.y, this.cam.dir.z], [this.cam.up.x, this.cam.up.y, this.cam.up.z]), [this.cam.dir.x, this.cam.dir.y, this.cam.dir.z]));
    this.cam.up.x = utmp[0]; this.cam.up.y = utmp[1]; this.cam.up.z = utmp[2];

    const ret = this.camRotKernel(
        this.move.tLR / 35,
        -this.move.tUD / 35,
        [this.cam.p.x, this.cam.p.y, this.cam.p.z],
        [this.cam.dir.x, this.cam.dir.y, this.cam.dir.z],
        [this.cam.up.x, this.cam.up.y, this.cam.up.z],
        1, 3.5, 2
    );

    this.cam.dir.x = ret[0];
    this.cam.dir.y = ret[1];
    this.cam.dir.z = ret[2];
    
    dlen = Math.sqrt(this.cam.dir.x*this.cam.dir.x + this.cam.dir.y*this.cam.dir.y + this.cam.dir.z*this.cam.dir.z);
    this.cam.p.x += (this.cam.dir.x/dlen) * this.move.tFR * dt * 6;
    this.cam.p.y += (this.cam.dir.y/dlen) * this.move.tFR * dt * 6;
    this.cam.p.z += (this.cam.dir.z/dlen) * this.move.tFR * dt * 6;

};