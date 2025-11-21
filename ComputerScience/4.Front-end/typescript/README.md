# TypeScript 介紹與基礎教學

[TypeScript](https://www.typescriptlang.org/) 是 JavaScript 的超集，由 Microsoft 開發，提供靜態型別與型別檢查功能，能提升開發效能與程式碼品質，特別適用於大型或多人協作的專案。

## 為什麼要使用 TypeScript？

* ✅ 靜態型別檢查，減少執行期錯誤
* ✅ 更佳的 IDE 提示與自動完成（IntelliSense）
* ✅ 更容易維護與重構
* ✅ 相容所有 JavaScript 程式碼

## 安裝與初始化

### 全域安裝 TypeScript 編譯器

```bash
npm install -g typescript
```

### 建立 tsconfig.json

```bash
tsc --init
```

這個設定檔可控制編譯行為，例如輸出目錄、是否嚴格型別檢查等。

## 基本語法範例

### 型別註記

```ts
let age: number = 25;
let name: string = "Alice";
let isActive: boolean = true;
```

### 陣列與物件

```ts
let scores: number[] = [80, 90, 100];

let user: { name: string; age: number } = {
  name: "Bob",
  age: 30
};
```

### 函式參數與回傳值型別

```ts
function add(x: number, y: number): number {
  return x + y;
}

function logMessage(msg: string): void {
  console.log(msg);
}
```

### 介面（Interface）與物件型別

```ts
interface Person {
  name: string;
  age: number;
  address?: string; // 可選屬性
}

const p: Person = { name: "Tom", age: 20 };
```

### 類別（Class）與修飾字

```ts
class Animal {
  // public 預設，外部可存取
  public name: string;

  // private 僅限類別內部存取
  private age: number;

  // protected 可在子類別存取
  protected type: string = "mammal";

  constructor(name: string, age: number) {
    this.name = name;
    this.age = age;
  }

  public speak(): void {
    console.log(`${this.name} makes a sound.`);
  }

  private getAge(): number {
    return this.age;
  }
}
```

## 編譯與執行

### 編譯 .ts 檔為 .js

```bash
tsc app.ts
```

### 即時執行 TypeScript（使用 ts-node）

```bash
npm install -g ts-node

ts-node app.ts
```

## 搭配 React、Vue 等框架

* **React**：使用 `.tsx` 檔案，可提供元件型別支援
* **Vue**：支援 `<script lang="ts">`，需搭配 `vue-tsc`

## 與 JavaScript 的差異

| 功能        | JavaScript | TypeScript |
| --------- | ---------- | ---------- |
| 靜態型別檢查    | ❌          | ✅          |
| 編譯階段錯誤    | ❌          | ✅          |
| 介面定義      | ❌          | ✅          |
| 支援 ES 新語法 | ✅          | ✅          |

---

TypeScript 是現代 JavaScript 開發的進階利器，適合追求穩定性、可維護性與大型專案協作的開發者。無論是前端或後端（如 Node.js），都值得學習並導入。
