# CSS 互動偽類（Pseudo-classes）介紹

CSS 偽類（pseudo-classes）可針對元素在特定狀態下套用樣式，常用於設計互動效果、使用者回饋與表單控制。

## 1. `:hover`

當滑鼠游標懸停在元素上時觸發。

```css
a:hover {
  color: red;
  text-decoration: underline;
}
```

## 2. `:active`

當元素正在被點擊（按下但未放開）時觸發。

```css
button:active {
  background-color: yellow;
  transform: scale(0.98);
}
```

## 3. `:focus`

當元素獲得焦點（如輸入框被點選）時觸發。

```css
input:focus {
  outline: none;
  border: 2px solid blue;
}
```

## 4. `:checked`

用於表單元素（checkbox 或 radio）被選取時。

```css
input:checked + label {
  font-weight: bold;
  color: green;
}
```

## 5. `:disabled` / `:enabled`

根據元素是否可互動來設定樣式。

```css
input:disabled {
  background-color: #eee;
  cursor: not-allowed;
}

input:enabled {
  background-color: white;
}
```

## 6. `:nth-child()` / `:first-child` / `:last-child`

針對元素在父元素中的相對位置。

```css
li:nth-child(odd) {
  background-color: #f9f9f9;
}

li:first-child {
  font-weight: bold;
}
```

---

這些偽類可大幅提升網頁的互動性與可用性，是設計使用者介面時不可或缺的重要工具。
