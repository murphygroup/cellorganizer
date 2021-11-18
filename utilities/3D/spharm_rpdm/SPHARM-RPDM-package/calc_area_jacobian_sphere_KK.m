function J_area_sph = calc_area_jacobian_sphere_KK(x) 
% was:
%(x11,x12,x13,x14,x21,x22,x23,x24,x31,x32,x33,x34)
%CALC_AREA_JACOBIAN_SPHERE
%    J_AREA_SPH = CALC_AREA_JACOBIAN_SPHERE(X11,X12,X13,X14,X21,X22,X23,X24,X31,X32,X33,X34)
%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    10-Mar-2018 21:49:15
% ------------------------------------------------------------------------
% Modifications June 2021 by Khaled Khairy: khaled.khairy@stjude.org
% St. Jude Children's Research Hospital
% Comments in code: Search for "KK"
%
% KK: modified to make vectorized for all faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(x,3);
x11 = x(1,1,:);
x12 = x(1,2,:);
x13 = x(1,3,:);
x14 = x(1,4,:);

x21 = x(2,1,:);
x22 = x(2,2,:);
x23 = x(2,3,:);
x24 = x(2,4,:);

x31 = x(3,1,:);
x32 = x(3,2,:);
x33 = x(3,3,:);
x34 = x(3,4,:);


t12 = x11.*x22.*x34;
t13 = x11.*x24.*x32;
t14 = x12.*x21.*x34;
t15 = x12.*x24.*x31;
t16 = x14.*x21.*x32;
t17 = x14.*x22.*x31;
t2 = t12-t13-t14+t15+t16-t17;
t4 = x11.*x12;
t5 = x21.*x22;
t6 = x31.*x32;
t7 = t4+t5+t6;
t8 = x11.*x14;
t9 = x21.*x24;
t10 = x31.*x34;
t11 = t8+t9+t10;
t18 = x12.*x14;
t19 = x22.*x24;
t20 = x32.*x34;
t21 = t7.*t11;
t3 = t18+t19+t20-t21;
t22 = t2.^2;
t23 = x12.*x13;
t24 = x22.*x23;
t25 = x32.*x33;
t26 = t23+t24+t25;
t27 = x11.*x22.*x33;
t28 = x12.*x23.*x31;
t29 = x13.*x21.*x32;
t31 = x11.*x23.*x32;
t32 = x12.*x21.*x33;
t33 = x13.*x22.*x31;
t30 = t27+t28+t29-t31-t32-t33;
t34 = x11.*x13;
t35 = x21.*x23;
t36 = x31.*x33;
t58 = t7.*t26;
t37 = t34+t35+t36-t58;
t38 = t30.^2;
t39 = x13.*x14;
t40 = x23.*x24;
t41 = x33.*x34;
t42 = t39+t40+t41;
t43 = x11.*x23.*x34;
t44 = x13.*x24.*x31;
t45 = x14.*x21.*x33;
t47 = x11.*x24.*x33;
t48 = x13.*x21.*x34;
t49 = x14.*x23.*x31;
t46 = t43+t44+t45-t47-t48-t49;
t64 = t11.*t42;
t50 = t34+t35+t36-t64;
t51 = t46.^2;
t52 = t3.^2;
t53 = t22+t52;
t54 = 1.0./t53;
t55 = 1.0./t2;
t56 = 1.0./t2.^2;
t57 = 1.0./t30;
t59 = 1.0./t30.^2;
t60 = t37.^2;
t61 = t38+t60;
t62 = 1.0./t61;
t63 = 1.0./t46;
t65 = 1.0./t46.^2;
t66 = t50.^2;
t67 = t51+t66;
t68 = 1.0./t67;
t69 = t7.*x14;
t70 = t11.*x12;
t71 = t69+t70;
t72 = t55.*t71;
t73 = x22.*x34;
t127 = x24.*x32;
t74 = t73-t127;
t75 = t3.*t56.*t74;
t76 = t72+t75;
t183 = t26.*x12;
t77 = -t183+x13;
t78 = t57.*t77;
t79 = x22.*x33;
t155 = x23.*x32;
t80 = t79-t155;
t184 = t37.*t59.*t80;
t81 = t78-t184;
t82 = t38.*t62.*t81;
t185 = t42.*x14;
t83 = -t185+x13;
t84 = t63.*t83;
t85 = x23.*x34;
t89 = x24.*x33;
t86 = t85-t89;
t186 = t50.*t65.*t86;
t87 = t84-t186;
t88 = t51.*t68.*t87;
t90 = x12.*x23.*x34;
t91 = x13.*x24.*x32;
t92 = x14.*x22.*x33;
t94 = x12.*x24.*x33;
t95 = x13.*x22.*x34;
t96 = x14.*x23.*x32;
t93 = t90+t91+t92-t94-t95-t96;
t102 = t26.*t42;
t97 = t18+t19+t20-t102;
t98 = t93.^2;
t99 = 1.0./t93;
t100 = x13.*x34;
t201 = x14.*x33;
t101 = t100-t201;
t103 = 1.0./t93.^2;
t104 = t97.^2;
t105 = t98+t104;
t106 = 1.0./t105;
t107 = x13.*x24;
t221 = x14.*x23;
t108 = t107-t221;
t109 = t7.*x13;
t110 = t26.*x11;
t111 = t109+t110;
t112 = t57.*t111;
t113 = x21.*x33;
t153 = x23.*x31;
t114 = t113-t153;
t230 = t37.*t59.*t114;
t115 = t112-t230;
t232 = t11.*x11;
t116 = -t232+x14;
t117 = t55.*t116;
t118 = x21.*x34;
t128 = x24.*x31;
t119 = t118-t128;
t120 = t3.*t56.*t119;
t121 = t117+t120;
t122 = t22.*t54.*t121;
t233 = t42.*x13;
t123 = -t233+x14;
t124 = t99.*t123;
t234 = t86.*t97.*t103;
t125 = t124-t234;
t126 = t98.*t106.*t125;
t129 = x12.*x34;
t191 = x14.*x32;
t130 = t129-t191;
t131 = x11.*x34;
t244 = x14.*x31;
t132 = t131-t244;
t133 = x12.*x24;
t211 = x14.*x22;
t134 = t133-t211;
t135 = x11.*x24;
t262 = x14.*x21;
t136 = t135-t262;
t137 = t26.*x14;
t138 = t42.*x12;
t139 = t137+t138;
t140 = t99.*t139;
t276 = t74.*t97.*t103;
t141 = t140-t276;
t278 = t7.*x12;
t142 = -t278+x11;
t143 = t57.*t142;
t144 = x21.*x32;
t154 = x22.*x31;
t145 = t144-t154;
t279 = t37.*t59.*t145;
t146 = t143-t279;
t147 = t38.*t62.*t146;
t280 = t11.*x14;
t148 = -t280+x11;
t149 = t63.*t148;
t150 = t50.*t65.*t119;
t151 = t149+t150;
t152 = t51.*t68.*t151;
t156 = x11.*x33;
t239 = x13.*x31;
t157 = t156-t239;
t158 = x11.*x32;
t289 = x12.*x31;
t159 = t158-t289;
t160 = x12.*x33;
t195 = x13.*x32;
t161 = t160-t195;
t162 = x11.*x23;
t258 = x13.*x21;
t163 = t162-t258;
t164 = x11.*x22;
t306 = x12.*x21;
t165 = t164-t306;
t166 = x12.*x23;
t216 = x13.*x22;
t167 = t166-t216;
t168 = t11.*x13;
t169 = t42.*x11;
t170 = t168+t169;
t171 = t63.*t170;
t172 = t50.*t65.*t114;
t173 = t171+t172;
t321 = t7.*x11;
t174 = -t321+x12;
t175 = t55.*t174;
t322 = t3.*t56.*t145;
t176 = t175-t322;
t177 = t22.*t54.*t176;
t323 = t26.*x13;
t178 = -t323+x12;
t179 = t99.*t178;
t324 = t80.*t97.*t103;
t180 = t179-t324;
t181 = t98.*t106.*t180;
t182 = t22.*t54.*t76;
t187 = t7.*x24;
t188 = t11.*x22;
t189 = t187+t188;
t190 = t55.*t189;
t226 = t3.*t56.*t130;
t192 = t190-t226;
t228 = t26.*x22;
t193 = -t228+x23;
t194 = t57.*t193;
t196 = t37.*t59.*t161;
t197 = t194+t196;
t198 = t38.*t62.*t197;
t229 = t42.*x24;
t199 = -t229+x23;
t200 = t63.*t199;
t202 = t50.*t65.*t101;
t203 = t200+t202;
t204 = t51.*t68.*t203;
t227 = t22.*t54.*t192;
t205 = t198+t204-t227;
t206 = t205.*x21;
t207 = t7.*x34;
t208 = t11.*x32;
t209 = t207+t208;
t210 = t55.*t209;
t212 = t3.*t56.*t134;
t213 = t210+t212;
t365 = t26.*x32;
t214 = -t365+x33;
t215 = t57.*t214;
t366 = t37.*t59.*t167;
t217 = t215-t366;
t218 = t38.*t62.*t217;
t367 = t42.*x34;
t219 = -t367+x33;
t220 = t63.*t219;
t368 = t50.*t65.*t108;
t222 = t220-t368;
t223 = t51.*t68.*t222;
t364 = t22.*t54.*t213;
t224 = t218+t223-t364;
t225 = t224.*x31;
t231 = t38.*t62.*t115;
t235 = t7.*x23;
t236 = t26.*x21;
t237 = t235+t236;
t238 = t57.*t237;
t240 = t37.*t59.*t157;
t241 = t238+t240;
t273 = t11.*x21;
t242 = -t273+x24;
t243 = t55.*t242;
t274 = t3.*t56.*t132;
t245 = t243-t274;
t246 = t22.*t54.*t245;
t275 = t42.*x23;
t247 = -t275+x24;
t248 = t99.*t247;
t249 = t97.*t101.*t103;
t250 = t248+t249;
t251 = t98.*t106.*t250;
t272 = t38.*t62.*t241;
t252 = t246+t251-t272;
t253 = t252.*x22;
t254 = t7.*x33;
t255 = t26.*x31;
t256 = t254+t255;
t257 = t57.*t256;
t372 = t37.*t59.*t163;
t259 = t257-t372;
t374 = t11.*x31;
t260 = -t374+x34;
t261 = t55.*t260;
t263 = t3.*t56.*t136;
t264 = t261+t263;
t265 = t22.*t54.*t264;
t375 = t42.*x33;
t266 = -t375+x34;
t267 = t99.*t266;
t376 = t97.*t103.*t108;
t268 = t267-t376;
t269 = t98.*t106.*t268;
t373 = t38.*t62.*t259;
t270 = t265+t269-t373;
t271 = t270.*x32;
t277 = t98.*t106.*t141;
t281 = t26.*x24;
t282 = t42.*x22;
t283 = t281+t282;
t284 = t99.*t283;
t285 = t97.*t103.*t130;
t286 = t284+t285;
t317 = t7.*x22;
t287 = -t317+x21;
t288 = t57.*t287;
t290 = t37.*t59.*t159;
t291 = t288+t290;
t292 = t38.*t62.*t291;
t318 = t11.*x24;
t293 = -t318+x21;
t294 = t63.*t293;
t319 = t50.*t65.*t132;
t295 = t294-t319;
t296 = t51.*t68.*t295;
t316 = t98.*t106.*t286;
t297 = t292+t296-t316;
t298 = t297.*x23;
t299 = t26.*x34;
t300 = t42.*x32;
t301 = t299+t300;
t302 = t99.*t301;
t380 = t97.*t103.*t134;
t303 = t302-t380;
t382 = t7.*x32;
t304 = -t382+x31;
t305 = t57.*t304;
t383 = t37.*t59.*t165;
t307 = t305-t383;
t308 = t38.*t62.*t307;
t384 = t11.*x34;
t309 = -t384+x31;
t310 = t63.*t309;
t311 = t50.*t65.*t136;
t312 = t310+t311;
t313 = t51.*t68.*t312;
t381 = t98.*t106.*t303;
t314 = t308+t313-t381;
t315 = t314.*x33;
t320 = t51.*t68.*t173;
t325 = t11.*x23;
t326 = t42.*x21;
t327 = t325+t326;
t328 = t63.*t327;
t358 = t50.*t65.*t157;
t329 = t328-t358;
t360 = t7.*x21;
t330 = -t360+x22;
t331 = t55.*t330;
t332 = t3.*t56.*t159;
t333 = t331+t332;
t334 = t22.*t54.*t333;
t361 = t26.*x23;
t335 = -t361+x22;
t336 = t99.*t335;
t337 = t97.*t103.*t161;
t338 = t336+t337;
t339 = t98.*t106.*t338;
t359 = t51.*t68.*t329;
t340 = t334+t339-t359;
t341 = t340.*x24;
t342 = t11.*x33;
t343 = t42.*x31;
t344 = t342+t343;
t345 = t63.*t344;
t346 = t50.*t65.*t163;
t347 = t345+t346;
t389 = t7.*x31;
t348 = -t389+x32;
t349 = t55.*t348;
t390 = t3.*t56.*t165;
t350 = t349-t390;
t351 = t22.*t54.*t350;
t391 = t26.*x33;
t352 = -t391+x32;
t353 = t99.*t352;
t392 = t97.*t103.*t167;
t354 = t353-t392;
t355 = t98.*t106.*t354;
t388 = t51.*t68.*t347;
t356 = t351+t355-t388;
t357 = t356.*x34;
t362 = t82+t88-t182;
t363 = t362.*x11;
t369 = t206+t225+t363;
t370 = t122+t126-t231;
t371 = t370.*x12;
t377 = t253+t271+t371;
t378 = t147+t152-t277;
t379 = t378.*x13;
t385 = t298+t315+t379;
t386 = t177+t181-t320;
t387 = t386.*x14;
t393 = t341+t357+t387;
J_area_sph = reshape([-t82-t88+t182+x11.*(t206+t225+x11.*(t82+t88-t22.*t54.*t76)),-t198-t204+t227+t369.*x21,-t218-t223+t364+t369.*x31,-t122-t126+t231+x12.*(t253+t271+x12.*(t122+t126-t38.*t62.*t115)),-t246-t251+t272+t377.*x22,-t265-t269+t373+t377.*x32,-t147-t152+t277+x13.*(t298+t315+x13.*(t147+t152-t98.*t106.*t141)),-t292-t296+t316+t385.*x23,-t308-t313+t381+t385.*x33,-t177-t181+t320+x14.*(t341+t357+x14.*(t177+t181-t51.*t68.*t173)),-t334-t339+t359+t393.*x24,-t351-t355+t388+t393.*x34],[3,4, n]);