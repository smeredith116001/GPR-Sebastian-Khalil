//#define colorRGB vec3
#define colorRGBA vec4
//Color 
vec4 colorScale(in vec2 fragCoord, in vec2 resolution)
{
    //return vec4(0.1,0.2,0.3,0.3);
    // color is 0 to 1 and coord is 0 to res
    // R and B -> horiz
    // G and A -> vert
    //return vec4(fragCoord, 0.0,1.0);

    // R  horiz
    // G  vert
    // return vec4(fragCoord, 0.0,1.0);
    // final: red-green gradient
    vec2 uv = fragCoord / resolution;
    //return vec4(uv, 0.5 , 1.0);
    //return vec4(uv, 0.25 , 1.0);
    vec3 color = vec3(uv, 0.25);
    float alpha = 1.0;
    return colorRGBA(color, alpha);
}

const float lineScale = 40.0;
//Got inspiration/learning from this guys effect. 
//I took time to learn how he did/what he did and tried to make a simplifed version that I could actually make without just copy pasting his work 
//Still Credits to:TheOnlyaaa
//http://glslsandbox.com/e#40822.0
//Full Screen Effect
float Cube(vec2 pos)
{
    float radius = atan(pos.x, pos.y);
    float num = abs(pos.x) + abs(pos.y);
    //float num = length(pos);
    float curLine = floor(num * lineScale);
    float speed = sin(curLine * 24.3);
    float offset = fract(sin(speed));

    float c = step(.4, tan(radius + offset + speed * iTime));

    float rnd = offset;
    return c * rnd;
}

vec4 rotatingDiamonds(float a, float b, out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    vec4 finalOutput = vec4(Cube(uv), Cube(uv - a), Cube(uv), b);
    return finalOutput;
}

//Circle
vec4 circleEffect(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 r = 2.0 * vec2(fragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    float radius = 0.5;

    vec4  color = vec4(r, 0.1 + 0.3 * sin(iTime), 1.0);
    //vec4 pixel;
        //if((uv.x*uv.x)+(uv.y*uv.y) < (radius*radius))
        //{
        //  pixel = color;
    //fragColor = pixel;
    vec4 white = vec4(0.0);
    vec4 pixel;
    pixel = white;
    if (length(r) < radius * radius)
    {
        pixel = color;
    }


    fragColor = pixel;
    return vec4(fragColor);
}

//CheckerBoard
vec4 checkerBoard(out vec4 fragColor, in vec2 fragCoord)
{
    float sizeOf = 10.0;
    vec2 Pos = floor(fragCoord / sizeOf);
    float checker = mod(Pos.x + mod(Pos.y, 2.0), 2.0);
    fragColor = checker * vec4(1.0, 1.0, 1.0, 1.0);
    return fragColor;
}
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    //vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
        //fragColor = colorScale(fragCoord,iResolution.xy);
        //fragColor = checkerBoard(fragColor, fragCoord);
        //fragColor = circleEffect(fragColor, fragCoord);
    fragColor = rotatingDiamonds(0.003, 1.0, fragColor, fragCoord);
}


// R, G, B, A = [0, 1]
//fragColor = vec4(1.0,1.0,1.0,1.0);

// Normalized pixel coordinates (from 0 to 1)
//vec2 uv = fragCoord/iResolution.xy;

// Time varying pixel color
//vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));

// Output to screen
//fragColor = vec4(col,1.0);