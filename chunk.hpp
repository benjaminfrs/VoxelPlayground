#include "block.hpp"
#include <array>
#include <OpenGL/gl3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <noise/noise.h>
//#include <noise/noiseutils.h>
#include "noiseutils.cpp"
#include <vector>

class Chunk
{
public:
  Chunk(){}

  ~Chunk(){}

  static const int CHUNK_SIZE = 40;

  auto GetBlocks() { return m_blocks; };

  void renderBlocks(const Shader &shader, const glm::mat4 &model){
    for(int x{0}; x < CHUNK_SIZE; x++){
      for(int y{0}; y < CHUNK_SIZE; y++){
        for(int z{0}; z < CHUNK_SIZE; z++){
          if(m_blocks[x][y][z].isActive()){
            //glm::mat4 model{1.0f};
            //model = glm::translate(model, glm::vec3(x, y, z));
            //std::cout << x << std::endl;
            //shader.setMat4("model", model);
            shader.setMat4("model", glm::translate(model, glm::vec3(x, y, z)));
            glDrawArrays(GL_TRIANGLES, 0, 36);
          }
  }}}};

  void findSurroundedBlocks(){
    //std::array<std::array<std::array<int, CHUNK_SIZE>, CHUNK_SIZE>, CHUNK_SIZE> temp{{0}};
    std::vector<glm::vec3> t;
    for(int x{1}; x < CHUNK_SIZE-1; x++){
      for(int y{1}; y < CHUNK_SIZE-1; y++){
        for(int z{1}; z < CHUNK_SIZE-1; z++){
          if(m_blocks[x-1][y][z].isActive() &&
             m_blocks[x][y-1][z].isActive() &&
             m_blocks[x][y][z-1].isActive() &&
             m_blocks[x+1][y][z].isActive() &&
             m_blocks[x][y+1][z].isActive() &&
             m_blocks[x][y][z+1].isActive()){
            //std::cout << "Found surounded block: " << x << y << z << std::endl;
            t.push_back(glm::vec3(x, y, z));
          }
        }
      }
    }
    for(glm::vec3 v : t){
      m_blocks[v.x][v.y][v.z].SetActive(false);
    }
  }

  void setBlockVisibility(const utils::NoiseMap& height_map){
    for(int x = 0; x < CHUNK_SIZE; x++){
        for(int z = 0; z < CHUNK_SIZE; z++){
            m_blocks[x][0][z].SetActive(true);

            std::cout << height_map.GetValue(x, z) << std::endl;
            //float height = height_map.GetValue(x, z) * (CHUNK_SIZE - 1) + CHUNK_SIZE/2;
            float height = height_map.GetValue(x, z) * (CHUNK_SIZE - 1);
            std::cout << std::endl << height << std::endl << "-----";
            //std::cout << height << std::endl;
            for(int y=1; y < height && y < CHUNK_SIZE; y++){
              m_blocks[x][y][z].SetActive(true);
              //std::cout << height_map[x][z] << ": " << y << std::endl;
            }
        }
    }
    findSurroundedBlocks();
  }

  void debugMakeChunkVisible(){
    for(int x{0}; x < CHUNK_SIZE; x++){
      for(int y{0}; y < CHUNK_SIZE; y++){
        for(int z{0}; z < CHUNK_SIZE; z++){
          m_blocks[x][y][z].SetActive(true);
        }
      }
    }
  }


  friend std::ostream& operator<<(std::ostream& os, Chunk& chunk);

private:
  std::array<std::array<std::array<Block, CHUNK_SIZE>, CHUNK_SIZE>, CHUNK_SIZE> m_blocks; 
};

std::ostream& operator<<(std::ostream& os, Chunk& chunk)
{
  for(int x=0; x < chunk.CHUNK_SIZE; x++){
    for(int y=0; y < chunk.CHUNK_SIZE; y++){
      for(int z=0; z < chunk.CHUNK_SIZE; z++){
        os << "Block[" << x << "][" << y << "][" << z << "]: " << chunk.GetBlocks()[x][y][z].isActive() << std::endl;
  }}}
  return os;
}

class ChunkList
{
public:
  ChunkList(){
    for(int i=0;i<MAX_CHUNKS;i++){
      chunks.push_back(Chunk());
    }
    //myModule.SetFrequency(2);
    //myModule.SetPersistence(0.75);
    lowerXBound = 0.0;
    lowerZBound = 0.0;
    upperXBound = boundIncrement;
    upperZBound = boundIncrement;
    heightMapBuilder.SetBounds(lowerXBound, upperXBound, lowerZBound, upperZBound);
    heightMapBuilder.SetDestSize(chunks[0].CHUNK_SIZE, chunks[0].CHUNK_SIZE);
    heightMapBuilder.SetDestNoiseMap(heightMap);
    heightMapBuilder.SetSourceModule(myModule);
    heightMapBuilder.Build();
  }

  ~ChunkList(){}

  static const int MAX_CHUNKS = 5;
  static const int boundIncrement = 1;

  void renderChunks(const Shader &shader, const glm::mat4 &model){
    for(int i=0;i<MAX_CHUNKS;i++){
      chunks[i].renderBlocks(shader, glm::translate(model, glm::vec3(chunks[i].CHUNK_SIZE*i + (i*2), 0, 0)));
    }
  }

  void setChunkVisibility(const utils::NoiseMap& height_map){
    for(auto &ch : chunks){
      lowerXBound += boundIncrement;
      upperXBound += boundIncrement;
      heightMapBuilder.SetBounds(lowerXBound, upperXBound, lowerZBound, upperZBound);
      heightMapBuilder.Build();

      ch.setBlockVisibility(heightMap);
    }
  }

  void debugVis(){
    for(auto &ch : chunks){
      ch.debugMakeChunkVisible();
    }
  }
private:
  std::vector<Chunk> chunks;
  module::Perlin myModule;
  utils::NoiseMap heightMap;
  utils::NoiseMapBuilderPlane heightMapBuilder;
  double lowerXBound;
  double upperXBound;
  double lowerZBound;
  double upperZBound;
};
