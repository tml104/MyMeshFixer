// Copyright 2024 tml104
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

https://blog.iyatt.com/?p=5832
build后复制依赖库：ldd MyMeshFixer | awk '{print $3}' | xargs -i cp -L {} .
每次build后：需要复制MyMesh到pack文件夹里面，然后用patchelf命令将外部库的elf添加到MyMeshFixer里面
打包：zip -r pack.zip ./pack
    - -r：递归要打包的目录